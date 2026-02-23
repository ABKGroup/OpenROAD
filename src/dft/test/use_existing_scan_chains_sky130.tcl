source "helpers.tcl"

read_lef sky130hd/sky130hd.tlef
read_lef sky130hd/sky130_fd_sc_hd_merged.lef
read_liberty sky130hd/sky130_fd_sc_hd__tt_025C_1v80.lib

read_verilog scan_architect_sky130_with_scan_ports.v
link_design scan_architect

create_clock -name clock1 -period 2.0000 -waveform {0.0000 1.0000} [get_ports {clock1}]
create_clock -name clock2 -period 2.0000 -waveform {0.0000 1.0000} [get_ports {clock2}]

scan_replace

# Import a user-defined scan path via SCANDEF/DEF scan chains stored in ODB.
read_def -incremental use_existing_scan_chains_sky130.scandef

# Stitch using the imported ODB scan chains (no re-architecting).
set_dft_config -use_existing_scan_chains 1
execute_dft_plan

set verilog_file [make_result_file use_existing_scan_chains_sky130.v]
write_verilog $verilog_file
diff_files $verilog_file use_existing_scan_chains_sky130.vok

set def_file [make_result_file use_existing_scan_chains_sky130.def]
write_def $def_file
diff_files $def_file use_existing_scan_chains_sky130.defok
