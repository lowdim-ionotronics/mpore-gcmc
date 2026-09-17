# mpore-gcmc
GCMC simulation of slit and cylindrical metallic pores using mplib

## Introduction

The mpore-gcmc package runs GCMC simulations of ions in both slit and
cylindrical pore geometries. The ions are treated as hard spheres with
charge. The pore walls are modelled as homogeneous surfaces, and the
particle energies are calculated using the [mplib
library](https://github.com/lowdim-ionotronics/mplib).

The GCMC simulations include:
1. Translational moves
2. Widom insertion/deletion
3. Ion swap

each of whose relative probabilities can be specified during usage (see
the [Notes](#notes) on `--prob-trans`/`--prob-widom` below -- they aren't
what they look like at first glance).

## Installation

```
pip install -e .
```

installs the `mpore_gcmc` library (needs `numpy`, `scipy`, pulled in
automatically). This does **not** install `mplib_ctypes` -- that's a
separate prerequisite from the sibling
[mplib](https://github.com/lowdim-ionotronics/mplib) repo: build it
(`./configure && make && make install` there) and make sure its install
prefix's `lib/` is on `LD_LIBRARY_PATH` before running anything here.

## Source files

- `mpore_gcmc/__init__.py` -- the library: pore geometry (`pore`,
  `cylinder`, `slit`), Monte Carlo machinery (`mc_settings`,
  `mcfunctions`), simulation state (`state`), output helpers (`output`),
  restart pickling (`pickle_write`/`pickle_load`), and the pore-width
  consistency check (`resolve_pore_width`).
- `run_gcmc.py` -- CLI driver: one voltage, an explicit list, or a range
  (three `-u/--voltage` values are expanded into `start,stop,step` via
  `numpy.arange`). The ion configuration carries over between voltages
  within one invocation. Two independent resume mechanisms are available
  for long runs -- see [Notes](#notes).
- `analyze_gcmc.py` -- post-processes `run_gcmc.py`'s output into
  charge/capacitance/density/energy vs. voltage, and optional
  density profiles (see [Analysis](#analysis)).

Both CLI scripts parse arguments and hand off to the `mpore_gcmc` library
-- neither defines any physics of its own.

## Usage

Either pass parameters as command-line flags, or point `-C`/`--config` at
a JSON file with the same field names as the long option names below
(e.g. `{"pore_width_accessible": 6.0, "pore_type": "cyl", ...}`) -- JSON
values fill in any option not given on the command line, so a CLI flag
always overrides the same key in the config file. See
`examples/jcp2026/*.json` for real examples.

```
python run_gcmc.py -C my_config.json
```

Command-line arguments can also be read from a file, since the parser is
configured with `fromfile_prefix_chars='@'`:

```
python run_gcmc.py @input_arguments.txt
```

The file should contain command-line arguments in the same format as a
normal terminal invocation.

## Command-line options

- `-o`, `--output-prefix`: output file prefix
  - Default value: `"Case"`
  - Type: string

- `-d`, `--output-dir`: directory to write output files into (created if
  missing, via `chdir` -- so this also scopes `--auto-resume`'s
  `./cont.restart` auto-detection to this directory rather than the
  invocation's own working directory). `analyze_gcmc.py` takes the same
  option, reading/writing in the same directory.
  - Default value: not set (current directory)
  - Type: string (directory path)

- `-A`, `--pore-width-accessible`: accessible pore width (in Å) -- the
  **diameter** for `--pore-type cyl`, the **gap width** for `--pore-type
  slit`. This is the width available to ion centres, not the width to the
  wall-atom/carbon centres (see `--pore-width-nominal`/`--wall-atom-radius`
  below). This accessible-vs-nominal distinction is stated once in the
  paper, in the MD methods subsection (`pore_geom/main.tex:123`) rather
  than the GCMC one, as `\wpore = \wporei - 2\carbonr`.
  - Provide this and/or `--pore-width-nominal`; at least one is required.
  - Type: floating-point number

- `--pore-width-nominal`: nominal pore width (in Å), to the
  wall-atom/carbon centres (diameter for `cyl`, gap width for `slit`).
  Related to `--pore-width-accessible` by `accessible = nominal -
  2*wall_atom_radius` (for **both** pore types -- see `--wall-atom-radius`
  below for why the factor is the same `2*` for a slit's two walls and a
  cylinder's one wall). If both `--pore-width-accessible` and
  `--pore-width-nominal` are given, they must be consistent with
  `--wall-atom-radius` or the program exits with an error.
  - Type: floating-point number

- `-a`, `--ion-radii`: ion radii (in Å)
  - Required: yes
  - Type: one or more floating-point numbers
  - At least two ion radii must be supplied.

- `-w`, `--mu-bulk`: bulk chemical potential for each particle type (in
  eV) -- the free-energy cost of transferring an ion from the bulk
  reservoir into the pore, before any electrostatic (voltage) term is
  added.
  - Required: yes
  - Type: one or more floating-point numbers

- `-u`, `--voltage`: one or more voltages (in V). Required -- omitting it
  exits with an error.
  - If exactly three values are given, they're interpreted as `start
    stop step` and expanded into a range via `numpy.arange`.
  - If a different count of values is given, they're used literally, in
    order (e.g. `-u 0.1 0.3 0.7` simulates precisely those three
    voltages, not a range).
  - Type: one or more floating-point numbers

- `-n`, `--production-steps`: number of Monte Carlo production steps
  - Default value: `1000000`
  - Type: integer

- `-t`, `--thermalization-steps`: number of Monte Carlo thermalization steps
  - Default value: `500000`
  - Type: integer

- `-P`, `--pore-type`: pore type: `cyl` or `slit`
  - Required: yes
  - Type: string
  - Validated: any other value raises an error.

- `-R`, `--restart`: `1` to resume from a restart file matching
  `--output-prefix`; omit or `0` to start fresh. Mutually exclusive with
  `--auto-resume` (see [Notes](#notes) for the difference between the two).
  - Default value: not set (fresh start)
  - Type: integer (0 or 1)

- `--auto-resume`: auto-detect `./cont.restart` at startup and resume the
  voltage list from there, with fine-grained mid-thermalization/
  mid-production checkpointing -- for long unattended runs that may get
  killed and resubmitted. Mutually exclusive with `-R`/`--restart` (see
  [Notes](#notes)).
  - Default value: off
  - Type: flag (no value)

- `-srh`, `--skip-restart-head`: `1` to skip the
  first voltage in the list when restarting (the restart file read still
  corresponds to that first voltage).
  - Default value: not set
  - Type: integer (0 or 1)

- `-T`, `--temperature`: temperature in K
  - Required: yes
  - Type: floating-point number

- `-e`, `--eshift`: distance (in Å) to shift the electron centre from the
  pore-atom centre
  - Default value: `0.0`
  - Type: floating-point number

- `-p`, `--epsr`: permittivity of the medium
  - Default value: `2.5`
  - Type: floating-point number

- `-L`, `--pore-length`: periodic simulation length (in Å) -- the
  cylinder axial length for `--pore-type cyl`, the square lateral box
  side length for `--pore-type slit`.
  - Default value: `100.0`
  - Type: floating-point number

- `-W`, `--wall-atom-radius`: radius of wall atoms (in Å; e.g. the carbon
  radius for a CNT), relating `--pore-width-accessible` and
  `--pore-width-nominal`.
  - Default value: `0.0`
  - Type: floating-point number

- `-s`, `--stat-frequency`: frequency for writing statistics
  - Default value: `1000`
  - Type: floating-point number

- `-r`, `--restart-frequency`: frequency for writing restart files
  - Default value: `1000`
  - Type: floating-point number

- `-c`, `--coords-frequency`: frequency for dumping particle coordinates
  - Default value: `None` (coordinate dumping is disabled)
  - Type: floating-point number

- `-pt`, `--prob-trans`: Monte Carlo move-selection threshold for
  translation moves. See [Notes](#notes) -- this is a cumulative
  threshold, not a standalone probability.
  - Default value: `0.5`
  - Type: floating-point number

- `-pw`, `--prob-widom`: Monte Carlo move-selection threshold up to and
  including Widom insertion/deletion moves. See [Notes](#notes).
  - Default value: `1.0`
  - Type: floating-point number

- `-C`, `--config`: JSON file of parameters, using the same field names as
  the long option names above. Supplies a default for any option not
  given on the command line; explicit CLI flags always take precedence.
  - Type: string (file path)

## Example usage

```bash
python run_gcmc.py -C examples/jcp2026/cyl_wpore6.json
```

or, entirely via CLI flags:

```bash
python run_gcmc.py \
  --output-prefix Case01 \
  --pore-width-accessible 24.0 \
  --ion-radii 1.8 2.0 \
  --mu-bulk -0.25 -0.25 \
  --voltage 0.0 0.8 0.05 \
  --production-steps 1000000 \
  --thermalization-steps 500000 \
  --pore-type cyl \
  --temperature 300.0 \
  --pore-length 100.0
```

## Notes

- `--pore-width-accessible` (or `--pore-width-nominal`), `--ion-radii`,
  `--mu-bulk`, `--voltage`, `--pore-type`, and `--temperature` are
  effectively required. The program prints an error message and exits if
  they are not provided.
- `--voltage` must contain exactly three values (`start`, `stop`, `step`)
  to be treated as a range -- see the full behavior under
  `-u`/`--voltage` above.
- Supply ion radii and chemical potentials as space-separated values
  after their respective options.
- The first two ion species are assigned charges `-1` and `+1`. Any
  additional ion species are assigned charge `0`.
- **`--prob-trans`/`--prob-widom` are cumulative thresholds inside a
  single random draw** (`mc_step()` picks translation if
  `rand < prob_trans`, Widom insertion/deletion if `prob_trans <= rand <
  prob_widom`, and a swap move otherwise) -- not the three independent
  probabilities the option names suggest. To get standalone probabilities
  translation=0.5, Widom=0.4, swap=0.1 (as used for the paper's Fig. 2,
  `pore_geom/main.tex:179`), pass `--prob-trans 0.5 --prob-widom 0.9`
  (`0.5 + 0.4`), not `--prob-widom 0.4`.
- Two independent resume mechanisms, not to be combined:
  - Plain (neither flag given): restart is **not** automatic -- a fresh
    simulation always starts, even if a matching restart file exists.
  - `-R 1`: resume from a specific, manually-named
    `<output-prefix>_<voltage>.restart` on every voltage in the list
    (e.g. to deliberately re-run or extend one saved voltage).
  - `--auto-resume`: resume is automatic -- if `./cont.restart` exists in
    the run directory, the simulation resumes from it (independent of
    `--output-prefix`) and figures out where in the voltage list to pick
    up, with fine-grained mid-thermalization/mid-production
    checkpointing. A per-voltage snapshot `v_<voltage>.restart` is also
    written after each voltage completes.
  - In all cases, restart files are written periodically (every
    `--restart-frequency` steps) and once more at the end of each
    voltage's run.

## Output files

- `<prefix>_<voltage>.count`: appended every `--stat-frequency` steps --
  one row of per-species particle counts (a histogram over the ion
  species) per write.
- `<prefix>_<voltage>.coords`: appended every `--coords-frequency` steps
  (only if `-c`/`--coords-frequency` is given) -- each block starts with
  a `step:<i_step>` header line followed by `[ion_type, x, y, z]` rows for
  every ion.
- `<prefix>_<voltage>.restart` (`-R`) / `./cont.restart` and
  `v_<voltage>.restart` (`--auto-resume`): pickled simulation state, for
  resuming a run (see Notes above).

## Analysis

`analyze_gcmc.py` computes accumulated charge, differential capacitance,
ion density, and stored energy vs. voltage from a run's `.count` files,
plus an optional z-density (slit) or radial-density (cyl) profile from
`.coords` files. It reads the *same* `-C/--config` JSON used to launch
the run, so pore geometry and ion charges don't need to be re-supplied:

```
python analyze_gcmc.py -C my_config.json --profile
```

writes `<prefix>_analysis.dat` (one row per voltage found for that
prefix, columns documented in the file's own header along with units and
the resolved parameter set) and, with `--profile`,
`<prefix>_<voltage>_profile.dat` per voltage with a `.coords` file.
Units mirror `single-file`: charge in µC/cm², capacitance in µF/cm²,
using the accessible pore radius for the surface-area normalization
(cylinder) or the flat wall area (slit). The reported energy,
`E_pore(u) = ∫₀ᵘC(u')u'du'` (µJ/cm²), is the single pore/electrode
stored energy. For a *symmetric* EDLC (identical electrodes, symmetric
electrolyte), the full-cell energy at cell voltage `v=2u` is
`2×E_pore(u)` -- not computed here, since that assumption doesn't hold
for arbitrary (e.g. asymmetric `--ion-radii`) runs. See the script's
module docstring for a caveat about the charge/capacitance normalization
prefactors not yet being validated against a known reference.

## License

GPLv3 -- see [`LICENSE`](LICENSE).

## Acknowledgements

The development of this code was supported by the
[National Science Centre (NCN)](https://www.ncn.gov.pl/), Poland, under Grant No. 2021 430/40/Q/ST4/00160 and HPC resources provided by [Institute for Computational Physics (ICP), University of Stuttgart](https://www.icp.uni-stuttgart.de/)
