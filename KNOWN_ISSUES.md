# Known issues

Notes from preparing this repo for public release, kept for transparency.

## Fixed: confinement used the wrong radius when `--wall-atom-radius != 0`

Before the `mpore_gcmc_serial.py`/`mpore_gcmc_serial_c.py` duplication was
merged into `mpore_gcmc/__init__.py`, the two files disagreed on which
radius bounds ion motion. `cylinder`/`slit`'s `volume_free`, `box_bound`,
and `random_position` (and `state.__init__`'s `r_xy`) used `self.Rc` (the
wall-atom/carbon-centre radius) in `mpore_gcmc_serial.py`, but `self.Rex`
(the accessible radius) in `mpore_gcmc_serial_c.py`. Only `Rex` is
correct: the paper (`pore_geom/main.tex:123`) defines the *accessible*
pore width as the one ion confinement should use, with the carbon-centre
width reserved for the electrostatic/image-charge plane (already
correctly implemented via `Rel`, derived from `Rc`, in `u1`/`u2`).

The merged library now uses `Rex` everywhere for confinement, matching
`mpore_gcmc_serial_c.py`'s (already-correct) behavior. This bug was dormant
in every run made with `--wall-atom-radius 0.0` (the only value used in
`run.sh` and, per available records, in producing the paper's Fig. 2),
since `Rc == Rex` when `wall_atom_radius == 0`.

## Caveat: `C_coeff`/`q_coeff` use the wrong radius, but are unused

`cylinder`/`slit.__init__` also compute `C_coeff`/`q_coeff` (intended for
converting simulated charge into physical charge/capacitance units) using
`Rc` rather than `Rex` -- structurally the same accessible-vs-nominal
question as the bug above, and as a similar bug fixed in the sibling
`single-file` repo's charge/capacitance normalization.

This has **not** been changed, because it turns out not to matter for any
published result: `grep -rn "C_coeff\|q_coeff"` across the pre-merge
source confirms both are dead attributes -- assigned once in
`cylinder.__init__`/`slit.__init__` and never read anywhere else in
either file, nor in either head script. The paper's actual
charge/capacitance/energy analysis (`pore_geom/figures-src/1d_vs_2d/`)
was done in separate, external scripts operating on the raw `.count`
output, not via these fields. Flagged here for anyone who starts relying
on `C_coeff`/`q_coeff` going forward -- the correct convention (`Rc` vs
`Rex`) for that use case hasn't been confirmed, so don't assume either is
right without checking.

## Where the accessible-width/carbon-radius convention is defined

`pore_geom/main.tex`'s dedicated GCMC methods subsection (`Monte Carlo
simulations for 1D and 2D pores`, lines 170-179) does not itself define
the accessible-vs-nominal pore-width relationship. It's stated once,
generically, in the MD methods subsection (line 123): "This definition of
pore width is consistent with that used in the GCMC simulations." If
you're looking for the actual $\wpore = \wporei - 2\carbonr$ relation,
that's where it is.
