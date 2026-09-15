# mpore-gcmc
GCMC simulation of slit and cylindrical metallic pores using Python

## Introduction

The mpore-gcmc package runs gcmc simulations of ions in both slit and cylindrical pore geometries. The ions are treated as hard-spheres with charge. The pore walls are modelled as homogenous surfaces and the particle energies are calculated using the [mplib library](https://github.com/lowdim-ionotronics/mplib). 

The gcmc simulations include,
1. Translational moves 
2. Widom insertion/deletion 
3. Ion swap 

each of whose relative probabilities can be specified during usage.

## Source Files

The code is placed in the ```./code``` folder. 

The files with ```"head"``` in their filename are the top-level files which accepts the simulation parameters and initializes the base classes and sets up the simulation. The files without ```"head"``` in their filename include all the base classes used for the simulation. 

The files with ```"_c"``` in their filename (mpore_gcmc_head_serial_c.py + mpore_gcmc_serial_c.py) can be used to run simulations for a range of electrode potentials with a specified step size (eg: 0.0 to 8.0 with a step size of 0.5). Here the final gcmc configuration from the previous potential step is used as the initial configuration for the next step. The files without ```"_c"``` in their filename can only be used to run the simulations at a single electrode potential. 

## Dependencies

- [mplib library](https://github.com/lowdim-ionotronics/mplib)
- Generic Python libraries, specifically [Numpy](https://numpy.org/) and [Scipy](https://scipy.org/)


## Usage

Arguments can also be read from a file because the parser is configured with:

```python
fromfile_prefix_chars='@'
```

Prefix the argument-file name with `@`:

```bash
python simulation.py @input_arguments.txt
```

The file should contain command-line arguments in the same format as a normal terminal invocation.


## Command-line options

- `-o`, `--output-prefix`: output file prefix
  - Default value: `"Case"`
  - Type: string

- `-A`, `--tube-radius`: accessible tube radius (in Å)
  - Required: yes
  - Type: floating-point number

- `-a`, `--ion-radii`: ion radii (in Å)
  - Required: yes
  - Type: one or more floating-point numbers
  - At least two ion radii must be supplied.

- `-w`, `--transfer-energy`: transfer energy for each particle type (in eV)
  - Required: yes
  - Type: one or more floating-point numbers

- `-u`, `--voltage`: voltage range start, stop, and step (in V) for files with `_c`in their file name. Provide a single value otherwise.
  - Default value: `[0.0, 0.8, 0.05]` when the supplied range does not contain exactly three values when using files with `_c`in their filename. The other files throws an error in this option is unspecified. 
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

- `-T`, `--temperature`: temperature in K
  - Required: yes
  - Type: floating-point number

- `-e`, `--eshift`: distance (in Å) to shift the electron centre from the pore-atom centre
  - Default value: `0.0`
  - Type: floating-point number

- `-p`, `--epsr`: permittivity of the medium
  - Default value: `2.5`
  - Type: floating-point number

- `-L`, `--tube-length`: length of the tube; for a slit pore, this is the slit size
  - Default value: `100.0`
  - Type: floating-point number

- `-W`, `--wall-atom-radius`: radius of wall atoms to reduce from the pore width
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

- `-pt`, `--prob-trans`: probability for Monte Carlo translation moves
  - Default value: `0.5`
  - Type: floating-point number

- `-pw`, `--prob-widom`: probability for Monte Carlo Widom moves
  - Default value: `1.0`
  - Type: floating-point number

## Example Usage

```bash
python simulation.py \
  --output-prefix Case01 \
  --tube-radius 12.0 \
  --ion-radii 1.8 2.0 \
  --transfer-energy -0.25 -0.25 \
  --voltage 0.0 0.8 0.05 \
  --production-steps 1000000 \
  --thermalization-steps 500000 \
  --pore-type cyl \
  --temperature 300.0 \
  --tube-length 100.0
```

## Notes

- `--tube-radius`, `--ion-radii`, `--transfer-energy`, `--voltage`, `--pore-type`, and `--temperature` are effectively required. The program prints an error message and exits if they are not provided.
- `--voltage` must contain exactly three values: `start`, `stop`, and `step`. If the provided array does not have three values, the script replaces it with `[0.0, 0.8, 0.05]` when you are using the code with `_c` in their file name.
- Supply ion radii and transfer energies as space-separated values after their respective options.
- The first two ion species are assigned charges `-1` and `+1`. Any additional ion species are assigned charge `0`.
- The code writes resart files containing configurations at each `--restart-frequency` step. Running the original command in the presence of `.restart` file inside the run folder will restart the simulations from their saved configuration.  

## Acknowledgements

The development of this code was supported by the
[National Science Centre (NCN)](https://www.ncn.gov.pl/), Poland, under Grant No. 2021 430/40/Q/ST4/00160 and HPC resources provided by [Institute for Computational Physics (ICP), University of Stuttgart](https://www.icp.uni-stuttgart.de/)