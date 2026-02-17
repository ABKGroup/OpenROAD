source "helpers.tcl"

read_lef sky130hd/sky130hd.tlef
read_lef sky130hd/sky130_fd_sc_hd_merged.lef
read_liberty sky130hd/sky130_fd_sc_hd__tt_025C_1v80.lib

read_verilog scan_architect_sky130.v
link_design scan_architect

create_clock -name clock1 -period 2.0000 -waveform {0.0000 1.0000} [get_ports {clock1}]
create_clock -name clock2 -period 2.0000 -waveform {0.0000 1.0000} [get_ports {clock2}]

set constraints_file [make_result_file group_cross_polarity_infeasible_sky130.constraints]
set fp [open $constraints_file w]
puts $fp "group g1 ff1_clk1_rising ff1_clk1_falling"
close $fp

set_dft_config -chain_count 2 -clock_mixing clock_mix -polarity_mode strict \
  -scan_order_constraints_file $constraints_file

scan_replace

set rc [catch { report_dft_plan -verbose } msg]
if { $rc == 0 } {
  error "Expected group_cross_polarity_infeasible_sky130 to fail due to cross-polarity grouping"
}
puts "Caught expected error: $msg"

