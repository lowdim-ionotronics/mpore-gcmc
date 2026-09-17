# JCP2026 examples

Six configs reproducing the GCMC curves of `fig:1D_vs_2D` (accessible pore
diameters 6/8/10 Å, cylindrical and slit) from our upcoming JCP2026 paper
(`pore_geom/main.tex`). Run any of them with, e.g.:

```
python run_gcmc.py -C examples/jcp2026/cyl_wpore6.json
```

Each config's `output_dir` sends its output (`.count`/`.restart`) into
`examples/jcp2026/output/`, created automatically. Analyze with, e.g.:

```
python analyze_gcmc.py -C examples/jcp2026/cyl_wpore6.json
```

Each JSON's `_source` field documents exactly which `main.tex` line
numbers each parameter came from, and the `prob_trans`/`prob_widom`
cumulative-threshold convention (see `run_gcmc.py --help` or the main
README).

**These use small illustrative step counts** (`thermalization_steps:
10000`, `production_steps: 100000`) so they run in seconds/minutes as a
usage demo -- they will **not** reproduce publication-quality statistics.
The paper used `1e6` equilibration and `1e7`-`1e8` production steps
(`main.tex:179`); edit those two fields (or override with `-t`/`-n` on
the command line) to reproduce the paper's actual runs.
