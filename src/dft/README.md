# DFT: Design for Testing

The Design for Testing module in OpenROAD (`dft`) is an implementation of Design For Testing.
New nets and logic are added to allow IC designs to be tested for errors in manufacturing.
Physical imperfections can cause hard failures and variability can cause timing errors.

A simple DFT insertion consist of the following parts:

* A scan_in pin where the test patterns are shifted in.
* A scan_out pin where the test patterns are read from.
* Scan cells that replace flops with registers that allow for testing.
* One or more scan chains (shift registers created from your scan cells).
* A scan_enable pin to allow your design to enter and leave the test mode.

## Commands

```{note}
- Parameters in square brackets `[-param param]` are optional.
- Parameters without square brackets `-param2 param2` are required.
```

### Set DFT Config 

The command `set_dft_config` sets the DFT configuration variables.

```tcl
	set_dft_config 
	    [-max_length <int>]
	    [-chain_count <int>]
	    [-max_chains <int>]
	    [-max_imbalance <float>]
	    [-clock_mixing <string>]
	    [-scan_order_metric <string>]
	    [-scan_order_solver <string>]
	    [-scanopt_rounds <int>]
	    [-scanopt_seed <int>]
	    [-scanopt_time_limit <float>]
	    [-scanopt_temp_control <0|1>]
	    [-scanopt_t_div <float>]
	    [-vertical_weight <float>]
	    [-timing_setup_weight <float>]
	    [-timing_hold_weight <float>]
	    [-timing_critical_slack <float>]
	    [-exclude_shift_registers <bool>]
	    [-prefer_qbar <bool>]
	    [-shift_register_min_length <int>]
	    [-scan_order_constraints_file <path>]
	    [-scan_enable_name_pattern <string>]
	    [-scan_in_name_pattern <string>]
	    [-scan_out_name_pattern <string>]
	    [-insert_lockup <bool>]
	    [-lockup_cell_rising <string>]
	    [-lockup_cell_falling <string>]
	    [-lockup_in_pin <string>]
	    [-lockup_out_pin <string>]
	    [-lockup_clock_pin_rising <string>]
	    [-lockup_clock_pin_falling <string>]
	    [-timing_buffer_cell <string>]
	    [-timing_buffer_in_pin <string>]
	    [-timing_buffer_out_pin <string>]
```

#### Options

| Switch Name | Description |
| ---- | ---- |
| `-max_length` | Hard maximum number of bits per scan chain. When set, infeasible constraints (e.g., groups larger than `max_length`) will error. |
| `-chain_count` | Exact total number of scan chains across the design. This takes priority over `max_chains`/`max_length` chain-count inference. In `no_mix`, this total will be distributed across clock domains as needed (at least one chain per clock). |
| `-max_chains` | Maximum total number of scan chains across the design. In `no_mix`, this must be at least the number of clock domains. |
| `-max_imbalance` | Maximum allowed chain length imbalance (percent). Constraint: `max(bits)/min(bits) <= 1 + max_imbalance/100`. Default is `30`. |
| `-clock_mixing` | How scan cells are partitioned into chains by clock. `no_mix` (default) does not mix different clock domains in a chain. `clock_mix` mixes clock domains (requires lockup insertion between domains). |
| `-polarity_mode` | How scan cells of different edge polarity are handled within a chain. `mid` (default) allows mixed polarity, stitching falling-edge cells before rising-edge cells in each chain. `strict` forbids mixing polarities within a chain, requiring separate chains when both polarities are present. |
| `-scan_order_metric` | Metric for ordering scan cells within each chain. `PLACEMENT` uses scan-pin Manhattan distance. `PIN_TO_NET` uses pin-to-net distance to global-route guides (or detailed routes when present), falling back to placement distance. |
| `-scan_order_solver` | Scan ordering solver. `HEURISTIC` is greedy + local cleanup. `SCANOPT` is an iterated local search (default). `UCLA_SCANOPT` uses the UCLA ScanOptpack-010411 reference implementation (fixed begin/end only; placement metric). |
| `-scanopt_rounds` | Iteration budget for `-scan_order_solver SCANOPT` (default `500000`). |
| `-scanopt_seed` | Random seed for `-scan_order_solver SCANOPT` (default `1`). |
| `-scanopt_time_limit` | Total time budget (seconds) for `SCANOPT` ordering across all scan chains. OpenROAD splits the budget across chains to keep runtime bounded as chain count increases. `0` means unlimited. |
| `-scanopt_temp_control` | Enable temperature control (optional uphill acceptance) for `SCANOPT` (`0`/`1`). |
| `-scanopt_t_div` | Temperature divisor for `SCANOPT` temperature control (larger reduces uphill acceptance). |
| `-vertical_weight` | Preferred wiring direction tuning. Values `>1` penalize vertical movement more than horizontal. Default `1.0`. |
| `-timing_setup_weight` | Optional timing-aware ordering penalty weight (setup). |
| `-timing_hold_weight` | Optional timing-aware ordering penalty weight (hold). |
| `-timing_critical_slack` | Slack threshold for timing-aware penalties (`0` = only penalize negative slack). |
| `-exclude_shift_registers` | Automatically detect simple functional shift-register chains (direct Q→D connections) and exclude them from `scan_replace` and scan planning (`0`/`1`). |
| `-prefer_qbar` | Prefer using the complemented output (`QN`/`Q_N`) as scan-out when the library does not tag a scan-out port (`0`/`1`). This can reduce added load on functional `Q` nets, at the cost of introducing inversion in the scan path. |
| `-shift_register_min_length` | Minimum chain length to classify as a shift register when `-exclude_shift_registers` is enabled (must be `>= 2`). Default is `4`. |
| `-scan_order_constraints_file` | Path to a scan ordering constraints file (see “Scan Ordering Constraints File” below). |
| `-scan_enable_name_pattern` | A format pattern with one or less set of braces (`{}`) to use to find or create scan enable drivers during scan chain stitching. The braces, if found, will be set to `0` as DFT architectures typically use a single shift-enable for all scan chains. If an un-escaped forward slash (`/`) is found, instead of searching for and/or creating a top-level port, an instance's pin will be searched for instead where the part of the string preceding the `/` is interpreted as the instance name and part succeeding it will be interpreted as the pin's name. |
| `-scan_in_name_pattern` | A format pattern with one or less braces (`{}`) to use to find or create scan in drivers during scan chain stitching. The braces will be replaced with the chain's ordinal number (starting at `0`). If an un-escaped forward slash (`/`) is found, instead of searching for and/or creating a top-level port, an instance's pin will be searched for instead where the part of the string preceding the `/` is interpreted as the instance name and part succeeding it will be interpreted as the pin's name. |
| `-scan_out_name_pattern` | A format pattern with one or less braces (`{}`) to use to find or create scan in loads during scan chain stitching. The braces will be replaced with the chain's ordinal number (starting at `0`). If an un-escaped forward slash (`/`) is found, instead of searching for and/or creating a top-level port, an instance's pin will be searched for instead where the part of the string preceding the `/` is interpreted as the instance name and part succeeding it will be interpreted as the pin's name. |
| `-insert_lockup` | Enable lockup insertion between adjacent scan cells when their clock domain differs (`0`/`1`). `clock_mix` implies lockup insertion. |
| `-lockup_cell_rising` | Library cell name used for a lockup latch on rising-edge domains. |
| `-lockup_cell_falling` | Library cell name used for a lockup latch on falling-edge domains. |
| `-lockup_in_pin` | Lockup cell input pin name. |
| `-lockup_out_pin` | Lockup cell output pin name. |
| `-lockup_clock_pin_rising` | Lockup clock pin name for rising-edge domains. |
| `-lockup_clock_pin_falling` | Lockup clock pin name for falling-edge domains. |
| `-timing_buffer_cell` | Optional buffer cell inserted after timing-critical scan-out sources (based on STA slack). |
| `-timing_buffer_in_pin` | Timing buffer input pin name. |
| `-timing_buffer_out_pin` | Timing buffer output pin name. |

### Scan Ordering Constraints File

The constraints file allows naming chains, defining per-chain begin/end points, and
expressing “ScanOpt-style” constraints such as grouping (contiguity), strict subpaths,
directed adjacencies, and partial order.

Notes:
- `#` starts a comment.
- Coordinates are specified in **DBU** (the same units used by DEF/ODB).
- A terminal can be either a top-level port name or an `inst/pin` reference. If an
  instance name contains a literal `/`, escape it as `\\/`.

Supported directives:
- `default_priority <0..127>`
- `group [<name>] [<priority>] <inst|group...>`: keep members together (supports hierarchical groups).
- `path [<name>] [<priority>] <inst...>`: strict order (no interpolation; forms a fixed subpath).
- `fixed_edge <from> <to>`: directed adjacency.
- `before <a> <b>`: partial order (enforces `<a>` before `<b>`; cycles error).
- `chain <name> [begin <x> <y>|<port|inst/pin>] [end <x> <y>|<port|inst/pin>]`
  - `begin/end` points are included in the ordering objective (begin→first, last→end).
  - If `begin/end` are terminals, they also override the scan-in/scan-out endpoint names for that chain.
  - If `begin/end` are points, the scan-in/scan-out ports (from name patterns) will have their pin locations set to those points during stitching.
- `chain_begin <name> <x> <y>|<port|inst/pin>` (alias: `chainbegin`)
- `chain_end <name> <x> <y>|<port|inst/pin>` (alias: `chainend`)
- `assign <chain> <inst|group...>`: force items into a named chain.
- `exclude <inst|group...>`: exclude instances from `scan_replace` and scan planning.
- `exclude_instance_pattern <glob...>`: exclude instances by name pattern.
- `exclude_master_pattern <glob...>`: exclude cell masters by name pattern.

### Report DFT Config

Prints the current DFT configuration to be used by `report_dft_plan` and
`execute_dft_plan`.

```tcl
report_dft_config
```

### Scan replace

Replaces flipflops with equivalent scan flipflops. This will generally be called before
placement, as it changes the area of cells.

```tcl
scan_replace
```

### Report DFT Plan

Prints a preview of the scan chains that will be stitched by `execute_dft_plan`. Use
this command to iterate and try different DFT configurations. This command does
not perform any modification to the design, and should be run after `scan_replace`
and global placement.

```tcl
report_dft_plan
    [-verbose]
```

#### Options

| Switch Name | Description |
| ---- | ---- |
| `-verbose` | Shows more information about each one of the scan chains that will be created. |

### Execute DFT Plan

Architect scan chains and connect them up in a way that minimises wirelength. As
a result, this should be run after placement, and after `scan_replace`.

```tcl
execute_dft_plan
```

### Write SCANDEF

Writes a DEF-style `SCANCHAINS` section (often called “SCANDEF”) describing the
stitched scan chains. This uses the scan chain objects stored in ODB by
`execute_dft_plan`.

```tcl
write_scandef -file <path>
```

### Scan Optimization

Reorders scan chains using the latest placement information (without re-running
`scan_replace`). This re-stitches scan connections using the current DFT config
(`set_dft_config`), including scan signal name patterns.


```tcl
scan_opt
```

## Example scripts

This example creates scan chains with a max length of 10 bits:

```
set_dft_config -max_length 10 -clock_mixing no_mix
report_dft_config
scan_replace
# Run global placement...
report_dft_plan -verbose
execute_dft_plan
```

## Regression tests

There are a set of regression tests in `./test`. For more information, refer to this [section](../../README.md#regression-tests).

Simply run the following script:

```shell
./test/regression
```


## Limitations

* Scan-chain optimization is heuristic and still evolving.
* Scan endpoints are controlled via the `-scan_*_name_pattern` options; OpenROAD will reuse matching existing ports/pins or create new ones when needed.
* There is currently no way to provide a full user-defined scan path (only partial constraints like `path`/`fixed_edge`).
* We can only work with one bit cells.

## License

BSD 3-Clause License. See [LICENSE](../../LICENSE) file.
