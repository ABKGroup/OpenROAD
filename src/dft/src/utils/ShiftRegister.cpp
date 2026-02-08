// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2023-2025, The OpenROAD Authors

#include "ShiftRegister.hh"

#include <unordered_map>
#include <utility>
#include <vector>

#include "ScanArchitectConfig.hh"
#include "Utils.hh"
#include "db_sta/dbNetwork.hh"
#include "db_sta/dbSta.hh"
#include "odb/db.h"
#include "sta/FuncExpr.hh"
#include "sta/Liberty.hh"
#include "sta/Sequential.hh"
#include "utl/Logger.h"

namespace dft::utils {

namespace {

struct SeqNode
{
  odb::dbInst* inst = nullptr;
  odb::dbNet* clk_net = nullptr;
  odb::dbNet* data_net = nullptr;
  odb::dbNet* out_net = nullptr;
  odb::dbNet* out_inv_net = nullptr;
};

bool isEligibleSequential(odb::dbInst* inst,
                          const ScanArchitectConfig& config,
                          sta::dbNetwork* db_network,
                          sta::LibertyCell*& liberty_cell)
{
  if (inst == nullptr) {
    return false;
  }
  if (inst->isDoNotTouch() || inst->isHierarchical()) {
    return false;
  }
  if (config.isInstanceExcluded(inst->getName(), inst->getMaster()->getName())) {
    return false;
  }

  sta::Cell* master_cell = db_network->dbToSta(inst->getMaster());
  liberty_cell = master_cell != nullptr ? db_network->libertyCell(master_cell)
                                        : nullptr;
  if (liberty_cell == nullptr) {
    return false;
  }
  if (!liberty_cell->hasSequentials() || liberty_cell->isClockGate()) {
    return false;
  }
  return true;
}

std::optional<SeqNode> makeSeqNode(odb::dbInst* inst,
                                  sta::dbSta* sta,
                                  const ScanArchitectConfig& config)
{
  sta::dbNetwork* db_network = sta->getDbNetwork();
  if (db_network == nullptr) {
    return std::nullopt;
  }

  sta::LibertyCell* liberty_cell = nullptr;
  if (!isEligibleSequential(inst, config, db_network, liberty_cell)) {
    return std::nullopt;
  }

  const sta::SequentialSeq& sequentials = liberty_cell->sequentials();
  if (sequentials.empty()) {
    return std::nullopt;
  }

  const sta::Sequential* seq = sequentials.front();
  if (!seq->isRegister()) {
    return std::nullopt;
  }

  // Only consider "simple" registers where the next-state is directly driven by
  // a single data input pin. This intentionally skips scan flops where next
  // state is a mux expression over functional + scan pins.
  sta::FuncExpr* data_expr = seq->data();
  if (data_expr == nullptr || data_expr->op() != sta::FuncExpr::op_port) {
    return std::nullopt;
  }

  sta::LibertyPort* data_port = data_expr->port();
  if (data_port == nullptr) {
    return std::nullopt;
  }

  odb::dbITerm* data_iterm = inst->findITerm(data_port->name());
  if (data_iterm == nullptr) {
    return std::nullopt;
  }

  odb::dbNet* data_net = data_iterm->getNet();
  if (data_net == nullptr) {
    return std::nullopt;
  }

  auto netFromPort = [&](sta::LibertyPort* port) -> odb::dbNet* {
    if (port == nullptr) {
      return nullptr;
    }
    odb::dbITerm* iterm = inst->findITerm(port->name());
    return iterm != nullptr ? iterm->getNet() : nullptr;
  };

  // Derive outputs. Some Liberty parsers may represent the sequential output as
  // an internal state; prefer common Q/QN ports when present.
  odb::dbNet* q_net = nullptr;
  odb::dbNet* qn_net = nullptr;
  if (liberty_cell != nullptr) {
    q_net = netFromPort(liberty_cell->findLibertyPort("Q"));
    qn_net = netFromPort(liberty_cell->findLibertyPort("QN"));
  }

  odb::dbNet* out_net = q_net != nullptr ? q_net : qn_net;
  odb::dbNet* out_inv_net = qn_net;

  if (out_net == nullptr) {
    out_net = netFromPort(seq->output());
  }
  if (out_inv_net == nullptr) {
    out_inv_net = netFromPort(seq->outputInv());
  }
  if (out_net == nullptr) {
    return std::nullopt;
  }
  if (out_inv_net == out_net) {
    out_inv_net = nullptr;
  }

  // Require a single clock net for shift-register detection.
  sta::FuncExpr* clk_expr = seq->clock();
  if (clk_expr == nullptr) {
    return std::nullopt;
  }
  sta::LibertyPort* clk_port = nullptr;
  if (clk_expr->op() == sta::FuncExpr::op_port) {
    clk_port = clk_expr->port();
  } else if (clk_expr->op() == sta::FuncExpr::op_not && clk_expr->left() != nullptr
             && clk_expr->left()->op() == sta::FuncExpr::op_port) {
    clk_port = clk_expr->left()->port();
  }
  if (clk_port == nullptr) {
    return std::nullopt;
  }
  odb::dbITerm* clk_iterm = inst->findITerm(clk_port->name());
  if (clk_iterm == nullptr) {
    return std::nullopt;
  }
  odb::dbNet* clk_net = clk_iterm->getNet();
  if (clk_net == nullptr) {
    return std::nullopt;
  }

  SeqNode node;
  node.inst = inst;
  node.clk_net = clk_net;
  node.data_net = data_net;
  node.out_net = out_net;
  node.out_inv_net = out_inv_net;
  return node;
}

void collectSeqNodes(odb::dbBlock* block,
                     sta::dbSta* sta,
                     const ScanArchitectConfig& config,
                     std::vector<SeqNode>& nodes)
{
  for (odb::dbInst* inst : block->getInsts()) {
    if (auto node = makeSeqNode(inst, sta, config)) {
      nodes.push_back(*node);
    }
  }
  for (odb::dbBlock* child : block->getChildren()) {
    collectSeqNodes(child, sta, config, nodes);
  }
}

int countITerms(odb::dbNet* net)
{
  int count = 0;
  for (odb::dbITerm* iterm : net->getITerms()) {
    (void) iterm;
    ++count;
  }
  return count;
}

}  // namespace

ShiftRegisterDetectResult DetectShiftRegisters(odb::dbDatabase* db,
                                               sta::dbSta* sta,
                                               const ScanArchitectConfig& config,
                                               utl::Logger* logger)
{
  ShiftRegisterDetectResult result;

  if (db == nullptr || sta == nullptr) {
    return result;
  }

  odb::dbChip* chip = db->getChip();
  if (chip == nullptr || chip->getBlock() == nullptr) {
    return result;
  }

  const int min_len = config.getShiftRegisterMinLength();
  if (min_len < 2) {
    return result;
  }

  std::vector<SeqNode> nodes;
  collectSeqNodes(chip->getBlock(), sta, config, nodes);
  if (nodes.empty()) {
    return result;
  }

  // Build net -> {drivers, sinks} maps for direct Q->D connections.
  std::unordered_map<odb::dbNet*, std::vector<const SeqNode*>> drivers_by_net;
  std::unordered_map<odb::dbNet*, std::vector<const SeqNode*>> sinks_by_net;
  drivers_by_net.reserve(nodes.size());
  sinks_by_net.reserve(nodes.size());

  for (const SeqNode& node : nodes) {
    sinks_by_net[node.data_net].push_back(&node);
    drivers_by_net[node.out_net].push_back(&node);
    if (node.out_inv_net != nullptr && node.out_inv_net != node.out_net) {
      drivers_by_net[node.out_inv_net].push_back(&node);
    }
  }

  // Build a directed graph where each node has at most one predecessor and one
  // successor (linear shift-register chains only).
  std::unordered_map<odb::dbInst*, odb::dbInst*> succ;
  std::unordered_map<odb::dbInst*, odb::dbInst*> pred;
  succ.reserve(nodes.size());
  pred.reserve(nodes.size());

  for (const auto& [net, sinks] : sinks_by_net) {
    if (sinks.size() != 1) {
      continue;
    }
    const auto it = drivers_by_net.find(net);
    if (it == drivers_by_net.end() || it->second.size() != 1) {
      continue;
    }

    // Require a simple two-instance net (driver + sink). This avoids excluding
    // generic pipelines where Q also fans out to other logic.
    if (countITerms(net) != 2) {
      continue;
    }

    const SeqNode* dst = sinks.front();
    const SeqNode* src = it->second.front();

    if (src->inst == nullptr || dst->inst == nullptr || src->inst == dst->inst) {
      continue;
    }
    if (src->clk_net == nullptr || dst->clk_net == nullptr
        || src->clk_net != dst->clk_net) {
      continue;
    }

    if (succ.find(src->inst) != succ.end()) {
      continue;
    }
    if (pred.find(dst->inst) != pred.end()) {
      continue;
    }
    succ[src->inst] = dst->inst;
    pred[dst->inst] = src->inst;
  }

  if (succ.empty()) {
    return result;
  }

  std::unordered_set<odb::dbInst*> visited;
  visited.reserve(succ.size());

  auto record_chain = [&](const std::vector<odb::dbInst*>& chain) {
    if (static_cast<int>(chain.size()) < min_len) {
      return;
    }
    ++result.chains;
    for (odb::dbInst* inst : chain) {
      if (inst != nullptr) {
        result.instances.insert(inst->getName());
      }
    }
  };

  // Walk all paths that have a unique start (no predecessor).
  for (const auto& [start, _] : succ) {
    if (pred.find(start) != pred.end()) {
      continue;
    }
    if (visited.find(start) != visited.end()) {
      continue;
    }

    std::vector<odb::dbInst*> chain;
    odb::dbInst* cur = start;
    while (cur != nullptr) {
      const bool inserted = visited.insert(cur).second;
      if (!inserted) {
        break;
      }
      chain.push_back(cur);
      const auto it = succ.find(cur);
      if (it == succ.end()) {
        break;
      }
      cur = it->second;
    }
    record_chain(chain);
  }

  // Walk any remaining cycles.
  for (const auto& [start, _] : succ) {
    if (visited.find(start) != visited.end()) {
      continue;
    }

    std::vector<odb::dbInst*> chain;
    odb::dbInst* cur = start;
    while (cur != nullptr) {
      const bool inserted = visited.insert(cur).second;
      if (!inserted) {
        break;
      }
      chain.push_back(cur);
      const auto it = succ.find(cur);
      if (it == succ.end()) {
        break;
      }
      cur = it->second;
    }
    record_chain(chain);
  }

  if (logger != nullptr && !result.instances.empty()) {
    logger->info(utl::DFT,
                 250,
                 "Auto-excluded {} shift-register instance(s) in {} chain(s) "
                 "(min_length={})",
                 result.instances.size(),
                 result.chains,
                 min_len);
  }

  return result;
}

}  // namespace dft::utils
