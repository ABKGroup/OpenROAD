source "helpers.tcl"

read_lef sky130hd/sky130hd.tlef
read_lef sky130hd/sky130_fd_sc_hd_merged.lef
read_liberty sky130hd/sky130_fd_sc_hd__tt_025C_1v80.lib

read_verilog scan_architect_sky130.v
link_design scan_architect

create_clock -name clock1 -period 2.0000 -waveform {0.0000 1.0000} [get_ports {clock1}]
create_clock -name clock2 -period 2.0000 -waveform {0.0000 1.0000} [get_ports {clock2}]

# Infeasible: With NoMix, at least 4 chains are required (one per clock/edge pair)
# when max_length is 5, but max_chains=3 caps the total chain count.
set_dft_config -max_length 5 -max_chains 3

scan_replace

set rc [catch { report_dft_plan -verbose } msg]
if { $rc == 0 } {
  error "Expected max_total_chain_count_sky130 to fail due to infeasible max_chains"
}
puts "Caught expected error: $msg"
