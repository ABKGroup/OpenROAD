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

sta::define_cmd_args "buffer_scan_enable" { -buffer_cell buffer_cell
                                            [-max_fanout max_fanout]
                                            [-max_levels max_levels] }
proc buffer_scan_enable { args } {
  sta::parse_key_args "buffer_scan_enable" args \
    keys {-buffer_cell -max_fanout -max_levels} \
    flags {}

  sta::check_argc_eq0 "buffer_scan_enable" $args

  if { [ord::get_db_block] == "NULL" } {
    utl::error DFT 295 "No design block found."
  }
  if { ![info exists keys(-buffer_cell)] } {
    utl::error DFT 296 "Missing required -buffer_cell argument."
  }

  set max_fanout 64
  if { [info exists keys(-max_fanout)] } {
    set max_fanout $keys(-max_fanout)
    sta::check_positive_integer "-max_fanout" $max_fanout
  }

  set max_levels 3
  if { [info exists keys(-max_levels)] } {
    set max_levels $keys(-max_levels)
    sta::check_positive_integer "-max_levels" $max_levels
  }

  dft::buffer_scan_enable $keys(-buffer_cell) $max_fanout $max_levels
}

sta::define_cmd_args "write_scandef" { -file file }
proc write_scandef { args } {
  sta::parse_key_args "write_scandef" args \
    keys {-file} \
    flags {}

  sta::check_argc_eq0 "write_scandef" $args

  if { [ord::get_db_block] == "NULL" } {
    utl::error DFT 269 "No design block found."
  }
  if { ![info exists keys(-file)] } {
    utl::error DFT 270 "Missing required -file argument."
  }

  dft::write_scandef $keys(-file)
}

sta::define_cmd_args "set_dft_config" { [-max_length max_length]
                                        [-chain_count chain_count]
                                        [-max_chains max_chains]
                                        [-max_imbalance max_imbalance]
                                        [-clock_mixing clock_mixing]
                                        [-polarity_mode polarity_mode]
                                        [-scan_order_metric scan_order_metric]
                                        [-scan_order_solver scan_order_solver]
                                        [-scanopt_rounds scanopt_rounds]
                                        [-scanopt_seed scanopt_seed]
                                        [-scanopt_time_limit scanopt_time_limit]
                                        [-scanopt_temp_control scanopt_temp_control]
                                        [-scanopt_t_div scanopt_t_div]
                                        [-vertical_weight vertical_weight]
                                        [-timing_setup_weight timing_setup_weight]
                                        [-timing_hold_weight timing_hold_weight]
                                        [-timing_critical_slack timing_critical_slack]
                                        [-exclude_shift_registers exclude_shift_registers]
                                        [-prefer_qbar prefer_qbar]
                                        [-shift_register_min_length shift_register_min_length]
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
      -max_imbalance
      -clock_mixing
      -polarity_mode
      -scan_order_metric
      -scan_order_solver
      -scanopt_rounds
      -scanopt_seed
      -scanopt_time_limit
      -scanopt_temp_control
      -scanopt_t_div
      -vertical_weight
      -timing_setup_weight
      -timing_hold_weight
      -timing_critical_slack
      -exclude_shift_registers
      -prefer_qbar
      -shift_register_min_length
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

  if { [info exists keys(-max_imbalance)] } {
    set v $keys(-max_imbalance)
    if { ![string is double -strict $v] } {
      utl::error DFT 102 "Expected a floating-point value for -max_imbalance"
    }
    if { $v < 0.0 } {
      utl::error DFT 103 "Expected a non-negative value for -max_imbalance"
    }
    dft::set_dft_config_max_imbalance $v
  }

  if { [info exists keys(-clock_mixing)] } {
    set clock_mixing $keys(-clock_mixing)
    dft::set_dft_config_clock_mixing $clock_mixing
  }

  if { [info exists keys(-polarity_mode)] } {
    set polarity_mode $keys(-polarity_mode)
    dft::set_dft_config_polarity_mode $polarity_mode
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

  if { [info exists keys(-scanopt_time_limit)] } {
    set s $keys(-scanopt_time_limit)
    if { ![string is double -strict $s] } {
      utl::error DFT 204 "Expected a floating-point value for -scanopt_time_limit"
    }
    if { $s < 0.0 } {
      utl::error DFT 205 "Expected a non-negative value for -scanopt_time_limit"
    }
    dft::set_dft_config_scanopt_time_limit $s
  }

  if { [info exists keys(-scanopt_temp_control)] } {
    set v $keys(-scanopt_temp_control)
    if { ![string is integer -strict $v] || ($v != 0 && $v != 1) } {
      utl::error DFT 217 "Expected 0 or 1 for -scanopt_temp_control"
    }
    dft::set_dft_config_scanopt_temp_control $v
  }

  if { [info exists keys(-scanopt_t_div)] } {
    set v $keys(-scanopt_t_div)
    if { ![string is double -strict $v] } {
      utl::error DFT 221 "Expected a floating-point value for -scanopt_t_div"
    }
    if { $v <= 0.0 } {
      utl::error DFT 222 "Expected a positive value for -scanopt_t_div"
    }
    dft::set_dft_config_scanopt_t_div $v
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

  if { [info exists keys(-exclude_shift_registers)] } {
    set v $keys(-exclude_shift_registers)
    if { ![string is boolean -strict $v] } {
      utl::error DFT 251 "-exclude_shift_registers must be a boolean (0/1/true/false)"
    }
    dft::set_dft_config_exclude_shift_registers [expr {$v ? 1 : 0}]
  }

  if { [info exists keys(-prefer_qbar)] } {
    set v $keys(-prefer_qbar)
    if { ![string is boolean -strict $v] } {
      utl::error DFT 253 "-prefer_qbar must be a boolean (0/1/true/false)"
    }
    dft::set_dft_config_prefer_qbar [expr {$v ? 1 : 0}]
  }

  if { [info exists keys(-shift_register_min_length)] } {
    set n $keys(-shift_register_min_length)
    sta::check_positive_integer "-shift_register_min_length" $n
    if { $n < 2 } {
      utl::error DFT 252 "Expected -shift_register_min_length >= 2"
    }
    dft::set_dft_config_shift_register_min_length $n
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
