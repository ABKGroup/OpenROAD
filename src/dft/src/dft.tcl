# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2023-2025, The OpenROAD Authors

sta::define_cmd_args "report_dft_plan" {[-verbose]}

proc report_dft_plan { args } {
  sta::parse_key_args "report_dft_plan" args \
    keys {} \
    flags {-verbose}

  sta::check_argc_eq0 "report_dft_plan" $args

  if { [ord::get_db_block] == "NULL" } {
    utl::error DFT 1 "No design block found."
  }

  set verbose [info exists flags(-verbose)]

  dft::report_dft_plan $verbose
}

sta::define_cmd_args "scan_replace" { }
proc scan_replace { args } {
  sta::parse_key_args "scan_replace" args \
    keys {} flags {}

  if { [ord::get_db_block] == "NULL" } {
    utl::error DFT 8 "No design block found."
  }
  dft::scan_replace
}

sta::define_cmd_args "execute_dft_plan" {}
proc execute_dft_plan { args } {
  sta::parse_key_args "execute_dft_plan" args \
    keys {} \
    flags {}

  if { [ord::get_db_block] == "NULL" } {
    utl::error DFT 9 "No design block found."
  }
  dft::execute_dft_plan
}

sta::define_cmd_args "set_dft_config" { [-max_length max_length]
                                        [-chain_count chain_count]
                                        [-max_chains max_chains]
                                        [-clock_mixing clock_mixing]
                                        [-scan_order_metric scan_order_metric]
                                        [-scan_order_solver scan_order_solver]
                                        [-scanopt_rounds scanopt_rounds]
                                        [-scanopt_seed scanopt_seed]
                                        [-vertical_weight vertical_weight]
                                        [-timing_setup_weight timing_setup_weight]
                                        [-timing_hold_weight timing_hold_weight]
                                        [-timing_critical_slack timing_critical_slack]
                                        [-scan_order_constraints_file scan_order_constraints_file]
                                        [-scan_enable_name_pattern scan_enable_name_pattern]
                                        [-scan_in_name_pattern scan_in_name_pattern]
                                        [-scan_out_name_pattern scan_out_name_pattern]
                                        [-insert_lockup insert_lockup]
                                        [-lockup_cell_rising lockup_cell_rising]
                                        [-lockup_cell_falling lockup_cell_falling]
                                        [-lockup_in_pin lockup_in_pin]
                                        [-lockup_out_pin lockup_out_pin]
                                        [-lockup_clock_pin_rising lockup_clock_pin_rising]
                                        [-lockup_clock_pin_falling lockup_clock_pin_falling]
                                        [-timing_buffer_cell timing_buffer_cell]
                                        [-timing_buffer_in_pin timing_buffer_in_pin]
                                        [-timing_buffer_out_pin timing_buffer_out_pin]
                                        }
proc set_dft_config { args } {
  sta::parse_key_args "set_dft_config" args \
    keys {
      -max_length
      -chain_count
      -max_chains
      -clock_mixing
      -scan_order_metric
      -scan_order_solver
      -scanopt_rounds
      -scanopt_seed
      -vertical_weight
      -timing_setup_weight
      -timing_hold_weight
      -timing_critical_slack
      -scan_order_constraints_file
      -scan_enable_name_pattern
      -scan_in_name_pattern
      -scan_out_name_pattern
      -insert_lockup
      -lockup_cell_rising
      -lockup_cell_falling
      -lockup_in_pin
      -lockup_out_pin
      -lockup_clock_pin_rising
      -lockup_clock_pin_falling
      -timing_buffer_cell
      -timing_buffer_in_pin
      -timing_buffer_out_pin
    } \
    flags {}

  sta::check_argc_eq0 "set_dft_config" $args

  if { [info exists keys(-max_length)] } {
    set max_length $keys(-max_length)
    sta::check_positive_integer "-max_length" $max_length
    dft::set_dft_config_max_length $max_length
  }

  if { [info exists keys(-chain_count)] } {
    set chain_count $keys(-chain_count)
    sta::check_positive_integer "-chain_count" $chain_count
    dft::set_dft_config_chain_count $chain_count
  }

  if { [info exists keys(-max_chains)] } {
    set max_chains $keys(-max_chains)
    sta::check_positive_integer "-max_chains" $max_chains
    dft::set_dft_config_max_chains $max_chains
  }

  if { [info exists keys(-clock_mixing)] } {
    set clock_mixing $keys(-clock_mixing)
    dft::set_dft_config_clock_mixing $clock_mixing
  }

  if { [info exists keys(-scan_order_metric)] } {
    set metric $keys(-scan_order_metric)
    dft::set_dft_config_scan_order_metric $metric
  }

  if { [info exists keys(-scan_order_solver)] } {
    set solver $keys(-scan_order_solver)
    dft::set_dft_config_scan_order_solver $solver
  }

  if { [info exists keys(-scanopt_rounds)] } {
    set rounds $keys(-scanopt_rounds)
    sta::check_positive_integer "-scanopt_rounds" $rounds
    dft::set_dft_config_scanopt_rounds $rounds
  }

  if { [info exists keys(-scanopt_seed)] } {
    set seed $keys(-scanopt_seed)
    sta::check_positive_integer "-scanopt_seed" $seed
    dft::set_dft_config_scanopt_seed $seed
  }

  if { [info exists keys(-vertical_weight)] } {
    set w $keys(-vertical_weight)
    if { ![string is double -strict $w] } {
      utl::error DFT 98 "Expected a floating-point value for -vertical_weight"
    }
    if { $w <= 0.0 } {
      utl::error DFT 99 "Expected a positive value for -vertical_weight"
    }
    dft::set_dft_config_vertical_weight $w
  }

  foreach {flag setter} {
    -timing_setup_weight dft::set_dft_config_timing_setup_weight
    -timing_hold_weight dft::set_dft_config_timing_hold_weight
    -timing_critical_slack dft::set_dft_config_timing_critical_slack
  } {
    if { [info exists keys($flag)] } {
      set v $keys($flag)
      if { ![string is double -strict $v] } {
        utl::error DFT 100 "Expected a floating-point value for $flag"
      }
      if { $v < 0.0 } {
        utl::error DFT 101 "Expected a non-negative value for $flag"
      }
      $setter $v
    }
  }

  if { [info exists keys(-scan_order_constraints_file)] } {
    set path $keys(-scan_order_constraints_file)
    dft::set_dft_config_scan_order_constraints_file $path
  }

  foreach {flag signal} {
    -scan_enable_name_pattern "scan_enable"
    -scan_in_name_pattern "scan_in"
    -scan_out_name_pattern "scan_out"
  } {
    if { [info exists keys($flag)] } {
      dft::set_dft_config_scan_signal_name_pattern $signal $keys($flag)
    }
  }

  if { [info exists keys(-insert_lockup)] } {
    set insert_lockup $keys(-insert_lockup)
    if { ![string is boolean -strict $insert_lockup] } {
      utl::error DFT 60 "-insert_lockup must be a boolean (0/1/true/false)"
    }
    dft::set_dft_config_insert_lockup [expr {$insert_lockup ? 1 : 0}]
  }

  if { [info exists keys(-lockup_cell_rising)] } {
    dft::set_dft_config_lockup_cell_rising $keys(-lockup_cell_rising)
  }
  if { [info exists keys(-lockup_cell_falling)] } {
    dft::set_dft_config_lockup_cell_falling $keys(-lockup_cell_falling)
  }
  if { [info exists keys(-lockup_in_pin)] } {
    dft::set_dft_config_lockup_in_pin $keys(-lockup_in_pin)
  }
  if { [info exists keys(-lockup_out_pin)] } {
    dft::set_dft_config_lockup_out_pin $keys(-lockup_out_pin)
  }
  if { [info exists keys(-lockup_clock_pin_rising)] } {
    dft::set_dft_config_lockup_clock_pin_rising $keys(-lockup_clock_pin_rising)
  }
  if { [info exists keys(-lockup_clock_pin_falling)] } {
    dft::set_dft_config_lockup_clock_pin_falling $keys(-lockup_clock_pin_falling)
  }

  if { [info exists keys(-timing_buffer_cell)] } {
    dft::set_dft_config_timing_buffer_cell $keys(-timing_buffer_cell)
  }
  if { [info exists keys(-timing_buffer_in_pin)] } {
    dft::set_dft_config_timing_buffer_in_pin $keys(-timing_buffer_in_pin)
  }
  if { [info exists keys(-timing_buffer_out_pin)] } {
    dft::set_dft_config_timing_buffer_out_pin $keys(-timing_buffer_out_pin)
  }
}

sta::define_cmd_args "report_dft_config" { }
proc report_dft_config { args } {
  sta::parse_key_args "report_dft_config" args keys {} flags {}
  dft::report_dft_config
}


sta::define_cmd_args "scan_opt" { }
proc scan_opt { args } {
  sta::parse_key_args "scan_opt" args \
    keys {} flags {}

  if { [ord::get_db_block] == "NULL" } {
    utl::error DFT 13 "No design block found."
  }
  dft::scan_opt
}
