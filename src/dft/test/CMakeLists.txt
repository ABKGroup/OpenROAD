or_integration_tests(
  "dft"
  TESTS
    max_chain_count_sky130
    max_total_chain_count_sky130
    exclude_shift_registers_nangate45
    prefer_qbar_nangate45
    one_cell_sky130
    config_flags_sky130
    buffer_scan_enable_sky130
    split_multibit_scan_cells_sky130
    place_sort_sky130
    ucla_scanopt_place_sort_sky130
    group_cross_polarity_infeasible_sky130
    scan_architect_clock_mix_sky130
    scan_architect_no_mix_sky130
    scan_architect_no_mix_nangate45
    scan_architect_register_bank_no_clock_mix_sky130
    scandef_core_sky130
    scandef_sky130
    sub_modules_sky130
    use_existing_scan_chains_sky130
    write_scandef_sky130
)


find_package(Boost CONFIG)
add_subdirectory(cpp)
