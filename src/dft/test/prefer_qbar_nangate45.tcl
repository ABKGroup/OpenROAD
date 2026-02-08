source "helpers.tcl"

read_lef Nangate45/Nangate45_tech.lef
read_lef Nangate45/Nangate45_stdcell.lef
read_liberty Nangate45/Nangate45_fast.lib

read_verilog scan_architect_no_mix_nangate45.v
link_design scan_architect_no_mix_nangate45

create_clock -name clk -period 2.0000 -waveform {0.0000 1.0000} [get_ports clk]

set_dft_config -max_length 10 -clock_mixing no_mix -prefer_qbar 1

scan_replace
execute_dft_plan

set result_dir [make_result_dir]
set cwd [pwd]
cd $result_dir
set out_scandef prefer_qbar_nangate45.scandef
write_scandef -file $out_scandef
cd $cwd

diff_files [file join $result_dir $out_scandef] prefer_qbar_nangate45.scandefok

