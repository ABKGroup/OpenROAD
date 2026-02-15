source "helpers.tcl"

read_lef sky130hd/sky130hd.tlef
read_lef sky130hd/sky130_fd_sc_hd_merged.lef
read_liberty sky130hd/sky130_fd_sc_hd__tt_025C_1v80.lib

read_def one_cell_sky130.def

create_clock -name main_clock -period 2.0000 -waveform {0.0000 1.0000} [get_ports {clock}]

set_dft_config -max_length 10 \
  -split_multibit_scan_cells 1 \
  -error_on_power_domain_crossings 1

report_dft_config

scan_replace
report_dft_plan
execute_dft_plan
