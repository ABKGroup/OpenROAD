source "helpers.tcl"

read_lef sky130hd/sky130hd.tlef
read_lef sky130hd/sky130_fd_sc_hd_merged.lef
read_liberty sky130hd/sky130_fd_sc_hd__tt_025C_1v80.lib

read_verilog place_sort_sky130.v
link_design place_sort

create_clock -name main_clock -period 2.0000 -waveform {0.0000 1.0000} [get_ports {clock}]

set_dft_config -max_length 10

scan_replace

proc place_inst { inst x y } {
  set db_inst [[ord::get_db_block] findInst $inst]
  $db_inst setLocation $x $y
  $db_inst setOrient R0
  $db_inst setPlacementStatus PLACED
}

place_inst ff1_clk1_rising 1000 2000
place_inst ff2_clk1_rising 5000 6000
place_inst ff3_clk1_rising 9500 9000
place_inst ff4_clk1_rising 9000 9000
place_inst ff5_clk1_rising 1000 1000
place_inst ff6_clk1_rising 7000 7000
place_inst ff7_clk1_rising 7500 7000
place_inst ff8_clk1_rising 4000 3000
place_inst ff9_clk1_rising 8000 8000
place_inst ff10_clk1_rising 3000 3000

execute_dft_plan

set inserted [buffer_scan_enable -buffer_cell sky130_fd_sc_hd__buf_1 -max_fanout 3 -max_levels 3]
puts "buffer_scan_enable inserted=$inserted"
check "buffer_scan_enable inserts buffers" {expr {$inserted > 0}} 1

