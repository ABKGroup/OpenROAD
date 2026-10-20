source "helpers.tcl"

read_lef sky130hd/sky130hd.tlef
read_lef sky130hd/sky130_fd_sc_hd_merged.lef
read_liberty sky130hd/sky130_fd_sc_hd__tt_025C_1v80.lib

read_verilog scan_inserted_design_sky130.v
link_design scan_inserted_design

# Load an existing SCANDEF into ODB, then export it back out via the DFT command.
read_def -incremental scan_inserted_design_sky130.scandef

set result_dir [make_result_dir]
set cwd [pwd]
cd $result_dir
set out_scandef write_scandef_sky130.scandef
write_scandef -file $out_scandef
cd $cwd

diff_files [file join $result_dir $out_scandef] write_scandef_sky130.scandefok
