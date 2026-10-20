source "helpers.tcl"

read_lef sky130hd/sky130hd.tlef
read_lef split_multibit_scan_cells_sky130.lef
read_liberty split_multibit_scan_cells_sky130.lib

read_def split_multibit_scan_cells_sky130.def

create_clock -name main_clock -period 2.0000 -waveform {0.0000 1.0000} [get_ports {clock}]

set_dft_config -max_length 10 -split_multibit_scan_cells 1 -scan_order_solver ILS

report_dft_plan -verbose
execute_dft_plan

set result_dir [make_result_dir]
set cwd [pwd]
cd $result_dir
set out_scandef split_multibit_scan_cells_sky130.scandef
write_scandef -file $out_scandef
cd $cwd

diff_files [file join $result_dir $out_scandef] split_multibit_scan_cells_sky130.scandefok
