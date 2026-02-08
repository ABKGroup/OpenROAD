#include <memory>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "ClockDomain.hh"
#include "Opt.hh"
#include "ScanArchitectConfig.hh"
#include "ScanCell.hh"
#include "gtest/gtest.h"
#include "odb/db.h"
#include "odb/geom.h"
#include "utl/Logger.h"

namespace dft::test {
namespace {

odb::dbBTerm* makeBTermAt(odb::dbBlock* block,
                          odb::dbTechLayer* layer,
                          const std::string& name,
                          int x,
                          int y)
{
  odb::dbNet* net = block->findNet(name.c_str());
  if (!net) {
    net = odb::dbNet::create(block, name.c_str());
  }
  odb::dbBTerm* term = odb::dbBTerm::create(net, name.c_str());
  odb::dbBPin* pin = odb::dbBPin::create(term);
  pin->setPlacementStatus(odb::dbPlacementStatus::PLACED);
  odb::dbBox::create(pin, layer, x, y, x, y);
  return term;
}

class PlacedScanCell final : public ScanCell
{
 public:
  PlacedScanCell(const std::string& name,
                 std::unique_ptr<ClockDomain> clock_domain,
                 utl::Logger* logger,
                 odb::dbBTerm* scan_in,
                 odb::dbBTerm* scan_out)
      : ScanCell(name, std::move(clock_domain), logger),
        scan_in_(scan_in),
        scan_out_(scan_out)
  {
  }

  uint64_t getBits() const override { return 1; }
  void connectScanEnable(const ScanDriver& /*pin*/) const override {}
  void connectScanIn(const ScanDriver& /*pin*/) const override {}
  void connectScanOut(const ScanLoad& /*pin*/) const override {}

  ScanLoad getScanEnable() const override
  {
    return ScanLoad(static_cast<odb::dbBTerm*>(nullptr));
  }
  ScanLoad getScanIn() const override { return ScanLoad(scan_in_); }
  ScanDriver getScanOut() const override { return ScanDriver(scan_out_); }

  odb::Point getOrigin() const override { return odb::Point(); }
  bool isPlaced() const override { return true; }

 private:
  odb::dbBTerm* scan_in_;
  odb::dbBTerm* scan_out_;
};

TEST(TestScanOrderOpt, UsesScanOutToScanInDistance)
{
  utl::Logger* logger = new utl::Logger();

  odb::dbDatabase* db = odb::dbDatabase::create();
  odb::dbTech* tech = odb::dbTech::create(db, "tech");
  odb::dbTechLayer* layer
      = odb::dbTechLayer::create(tech, "M1", odb::dbTechLayerType::ROUTING);
  odb::dbLib::create(db, "lib", tech, ',');
  odb::dbChip* chip = odb::dbChip::create(db, tech);
  odb::dbBlock* block = odb::dbBlock::create(chip, "block");

  // Construct three scan cells where A's scan_in is closest to B's scan_in, but
  // A's scan_out is closest to C's scan_in. A placement-based optimizer that
  // correctly models scan stitching should choose A->C, not A->B.
  odb::dbBTerm* a_si = makeBTermAt(block, layer, "a_si", 0, 0);
  odb::dbBTerm* a_so = makeBTermAt(block, layer, "a_so", 100, 0);
  odb::dbBTerm* b_si = makeBTermAt(block, layer, "b_si", 10, 0);
  odb::dbBTerm* b_so = makeBTermAt(block, layer, "b_so", 10, 0);
  odb::dbBTerm* c_si = makeBTermAt(block, layer, "c_si", 110, 0);
  odb::dbBTerm* c_so = makeBTermAt(block, layer, "c_so", 110, 0);

  ScanArchitectConfig config;
  config.setScanOrderMetric(ScanArchitectConfig::ScanOrderMetric::Placement);
  config.setScanOrderSolver(ScanArchitectConfig::ScanOrderSolver::Heuristic);

  std::vector<std::unique_ptr<ScanCell>> cells;
  cells.push_back(std::make_unique<PlacedScanCell>(
      "A",
      std::make_unique<ClockDomain>("clk", ClockEdge::Rising),
      logger,
      a_si,
      a_so));
  cells.push_back(std::make_unique<PlacedScanCell>(
      "B",
      std::make_unique<ClockDomain>("clk", ClockEdge::Rising),
      logger,
      b_si,
      b_so));
  cells.push_back(std::make_unique<PlacedScanCell>(
      "C",
      std::make_unique<ClockDomain>("clk", ClockEdge::Rising),
      logger,
      c_si,
      c_so));

  // Fix the chain start near A so the heuristic can't improve QoR by selecting a
  // different break point in an otherwise equivalent cycle.
  ScanArchitectConfig::ChainEndpoints endpoints;
  ScanArchitectConfig::ChainEndpoint begin;
  begin.type = ScanArchitectConfig::ChainEndpoint::Type::Point;
  begin.point.x = 0;
  begin.point.y = 0;
  endpoints.begin = begin;

  OptimizeScanWirelength(cells, config, logger, endpoints);

  ASSERT_EQ(cells.size(), 3);
  EXPECT_EQ(cells[0]->getName(), std::string_view("A"));
  EXPECT_EQ(cells[1]->getName(), std::string_view("C"));
  EXPECT_EQ(cells[2]->getName(), std::string_view("B"));
}

}  // namespace
}  // namespace dft::test
