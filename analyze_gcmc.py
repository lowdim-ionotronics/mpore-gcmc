# Computes accumulated charge, differential capacitance, ion density, and
# stored energy vs. voltage from run_gcmc.py output, plus an optional
# density profile. Originally written as four separate
# scripts by Vishnu Prasad Kurupath
# (https://vishnu-prasad-kurupath.github.io/); consolidated here into one.
#
# Units: charge in uC/cm^2, capacitance in uF/cm^2, density in mol/L,
# energy in uJ/cm^2. For a 1D (cylindrical) pore, the *accessible* pore
# width is used for the surface-area normalization, consistent with
# lowdim-ionotronics/single-file. (The exact charge/capacitance
# normalization prefactors are dimensionally consistent but haven't been
# validated against a known reference -- sanity-check a real result
# before trusting it for publication.)
#
# Note: the reported energy, E_pore(u) = integral_0^u C(u')u'du', is the
# single pore/electrode stored energy. For a *symmetric* EDLC (two
# identical electrodes, symmetric electrolyte), the full-cell energy at
# cell voltage v=2u is 2*E_pore(u) -- not computed here, since that
# equal-voltage-split assumption doesn't hold in general (e.g. for
# asymmetric --ion-radii), and this script has no way to know whether it
# applies to a given run.

import sys
import os
import re
import glob
import json
import argparse
import numpy as np
from scipy.integrate import cumulative_trapezoid
from scipy import constants as sc

import mpore_gcmc

E_CHARGE = sc.e      # C
K_B = sc.k           # J/K
N_AVOGADRO = sc.N_A  # 1/mol
ANGSTROM3_PER_LITRE = 1.0e27  # 1 L = 1e27 A^3


def load_config(path):
    with open(path) as f:
        return json.load(f)


def build_pore(cfg):
    """Construct the mpore_gcmc pore object (cylinder or slit) implied by
    a run config, giving access to the exact geometry (Rex, length,
    rho_factor) the simulation itself used -- without re-running the
    (potentially expensive) initial ion-placement loop that a full
    mpore_gcmc.state(...) would trigger.

    Returns (pore_obj, aion, d_ion_A, q_comp).
    """
    ptype = cfg['pore_type']
    if ptype not in ('cyl', 'slit'):
        raise ValueError(f"Unknown pore_type '{ptype}': must be 'cyl' or 'slit'")

    wall_atom_radius = cfg.get('wall_atom_radius', 0.0)
    pore_width_accessible, _ = mpore_gcmc.resolve_pore_width(
        cfg.get('pore_width_accessible'), cfg.get('pore_width_nominal'),
        wall_atom_radius)

    aion = np.array(cfg['ion_radii'])
    d_ion_A = 2.0 * np.max(aion)
    q_comp = np.array([-1, 1])
    for _ in range(len(aion) - 2):
        q_comp = np.append(q_comp, 0)

    pore_cls = mpore_gcmc.cylinder if ptype == 'cyl' else mpore_gcmc.slit
    pore_obj = pore_cls(
        pore_width_accessible, cfg.get('pore_length', 100.0),
        cfg.get('eshift', 0.0), wall_atom_radius,
        aion, K_B, d_ion_A, E_CHARGE, cfg['temperature'])
    return pore_obj, aion, d_ion_A, q_comp


def areal_factors(pore_obj, ptype, d_ion_A, T):
    """(q_factor, cap_factor): prefactors converting a dimensionless
    charge count (or charge-fluctuation variance) into uC/cm^2 (uF/cm^2),
    for use as q_factor * pore_obj.rho_factor * <raw quantity>. See the
    module docstring's CAVEAT.
    """
    d_ion_m = d_ion_A * 1e-10
    if ptype == 'cyl':
        Rex_m = pore_obj.Rex * 1e-10
        q_factor = E_CHARGE / (2.0 * np.pi * Rex_m * d_ion_m) * 100.0
        cap_factor = E_CHARGE**2 / (2.0 * np.pi * Rex_m * d_ion_m * K_B * T) * 100.0
    else:
        q_factor = E_CHARGE / (d_ion_m**2) * 100.0
        cap_factor = E_CHARGE**2 / (d_ion_m**2 * K_B * T) * 100.0
    return q_factor, cap_factor


def pore_volume_A3(pore_obj, ptype):
    if ptype == 'cyl':
        return np.pi * pore_obj.Rex**2 * pore_obj.length
    return pore_obj.length_x * pore_obj.length_y * (2.0 * pore_obj.Rex)


def find_count_files(prefix):
    """Return [(voltage, path), ...] sorted by voltage, for every
    <prefix>_<voltage>.count file present in the current directory
    (prefix may itself include a directory component, e.g. 'out/run')."""
    pattern = f"{prefix}_*.count"
    prefix_basename = os.path.basename(prefix)
    found = []
    for path in glob.glob(pattern):
        m = re.fullmatch(re.escape(prefix_basename) + r"_(.+)\.count", os.path.basename(path))
        if m:
            try:
                voltage = float(m.group(1))
            except ValueError:
                continue
            found.append((voltage, path))
    found.sort(key=lambda t: t[0])
    return found


def load_count_data(path):
    data = np.loadtxt(path)
    return np.atleast_2d(data)


def charge_and_capacitance(count_data, q_comp):
    """(mean_nq, var_nq): dimensionless mean charge and charge-fluctuation
    variance <(sum n_i q_i)^2> - <sum n_i q_i>^2, from raw per-species MC
    ion counts (rows = MC snapshots, columns = species)."""
    mean_counts = np.mean(count_data, axis=0)
    mean_nq = np.sum(mean_counts * q_comp)
    correl = np.mean(count_data[:, :, np.newaxis] * count_data[:, np.newaxis, :], axis=0)
    q_correl = np.outer(q_comp, q_comp)
    nq2_mean = np.sum(correl * q_correl)
    var_nq = nq2_mean - mean_nq**2
    return mean_nq, var_nq


def read_coords_frames(path):
    """Parse a .coords file written by mpore_gcmc.output.dump_coords:
    blocks starting with a 'step:<n>' header line, followed by
    [ion_type, x, y, z] rows. Returns a list of (n_ions, 4) arrays."""
    frames = []
    current = []
    with open(path) as f:
        for line in f:
            if line.startswith('step:'):
                if current:
                    frames.append(np.array(current, dtype=float))
                current = []
            else:
                current.append([float(x) for x in line.split()])
    if current:
        frames.append(np.array(current, dtype=float))
    return frames


def compute_profile(pore_obj, ptype, frames, nbins):
    """z-density profile (slit: confined z-direction, bins over
    +/-Rex) or radial density profile (cyl: bins over [0, Rex])."""
    if ptype == 'slit':
        edges = np.linspace(-pore_obj.Rex, pore_obj.Rex, nbins + 1)
        coord_index = 3  # ion_type, x, y, z
    else:
        edges = np.linspace(0.0, pore_obj.Rex, nbins + 1)
        coord_index = None  # computed below as sqrt(x^2+y^2)
    centers = 0.5 * (edges[1:] + edges[:-1])
    hist_sum = np.zeros(nbins)
    for frame in frames:
        if frame.size == 0:
            continue
        if ptype == 'slit':
            values = frame[:, coord_index]
        else:
            values = np.sqrt(frame[:, 1]**2 + frame[:, 2]**2)
        h, _ = np.histogram(values, bins=edges)
        hist_sum += h
    mean_hist = hist_sum / max(len(frames), 1)
    return centers, mean_hist


def main():
    parser = argparse.ArgumentParser(
        description="Analyze run_gcmc.py output: charge, capacitance, ion "
                     "density, and stored energy vs. voltage (from .count "
                     "files), plus an optional z/radial density profile "
                     "(from .coords files).")
    parser.add_argument("-C", "--config", required=True, metavar="FILE",
                         help="the run's own JSON config (same file passed to "
                              "run_gcmc.py via -C)")
    parser.add_argument("-o", "--output-prefix", metavar="FILE",
                         help="output-prefix to read <prefix>_<voltage>.count "
                              "(and, with --profile, .coords) files for; "
                              "defaults to the config's own output_prefix")
    parser.add_argument("-d", "--output-dir", metavar="DIR",
                         help="directory the run's output files live in (and "
                              "where this script's own output is written); "
                              "defaults to the config's own output_dir, else "
                              "the current directory")
    parser.add_argument("--out", metavar="FILE",
                         help="analysis output file; defaults to "
                              "<prefix>_analysis.dat")
    parser.add_argument("--profile", action="store_true",
                         help="also compute a z-density (slit) or radial "
                              "density (cyl) profile from <prefix>_<voltage>.coords "
                              "files, written to <prefix>_<voltage>_profile.dat")
    parser.add_argument("--profile-bins", type=int, default=100, metavar="N",
                         help="number of bins for --profile (default: 100)")
    args = parser.parse_args()

    cfg = load_config(args.config)

    output_dir = args.output_dir or cfg.get('output_dir')
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
        os.chdir(output_dir)

    prefix = args.output_prefix or cfg.get('output_prefix', 'Case')
    prefix_basename = os.path.basename(prefix)
    ptype = cfg['pore_type']

    pore_obj, aion, d_ion_A, q_comp = build_pore(cfg)
    q_factor, cap_factor = areal_factors(pore_obj, ptype, d_ion_A, cfg['temperature'])
    volume_A3 = pore_volume_A3(pore_obj, ptype)

    count_files = find_count_files(prefix)
    if not count_files:
        print(f"No {prefix}_*.count files found in the current directory.")
        sys.exit(1)

    u_list, Q_list, C_list, dens_list = [], [], [], []
    for u, path in count_files:
        count_data = load_count_data(path)
        mean_nq, var_nq = charge_and_capacitance(count_data, q_comp)
        Q_uC_cm2 = q_factor * pore_obj.rho_factor * mean_nq
        C_uF_cm2 = cap_factor * pore_obj.rho_factor * var_nq
        mean_total = np.sum(np.mean(count_data, axis=0))
        dens_mol_L = mean_total / (volume_A3 / ANGSTROM3_PER_LITRE) / N_AVOGADRO
        u_list.append(u)
        Q_list.append(Q_uC_cm2)
        C_list.append(C_uF_cm2)
        dens_list.append(dens_mol_L)

    u_arr = np.array(u_list)
    Q_arr = np.array(Q_list)
    C_arr = np.array(C_list)
    dens_arr = np.array(dens_list)

    # Single pore/electrode stored energy, E_pore(u) = int_0^u C(u')u'du'
    # -- see the module docstring's Note on the (not computed here)
    # symmetric-EDLC full-cell energy.
    E_pore_arr = cumulative_trapezoid(C_arr * u_arr, u_arr, initial=0.0)

    outname = args.out or f"{prefix}_analysis.dat"
    with open(outname, 'w') as f:
        f.write(f"# Analysis of {prefix}_<voltage>.count via analyze_gcmc.py\n")
        f.write(f"# Source files: {', '.join(p for _, p in count_files)}\n")
        f.write(f"# Config: {os.path.abspath(args.config)}\n")
        f.write(f"# pore_type={ptype}  T={cfg['temperature']}K  epsr={cfg['epsr']}  "
                f"pore_width_accessible={pore_obj.Rex*2:.4f}A  "
                f"wall_atom_radius={cfg.get('wall_atom_radius', 0.0)}A  "
                f"pore_length={cfg.get('pore_length', 100.0)}A  "
                f"ion_radii={[float(a) for a in aion]}\n")
        f.write("# u = electrode potential [V]\n")
        f.write("# E_pore = single pore/electrode energy = int C(u')u'du' [uJ/cm2]. "
                "For a symmetric EDLC (identical electrodes, symmetric "
                "electrolyte), the full-cell energy at cell voltage v=2u is "
                "2*E_pore(u).\n")
        f.write("# (1)u[V]  (2)Q[uC/cm2]  (3)C[uF/cm2]  (4)dens[mol/L]  "
                "(5)E_pore[uJ/cm2]\n")
        for row in zip(u_arr, Q_arr, C_arr, dens_arr, E_pore_arr):
            f.write("  ".join(f"{x:.6e}" for x in row) + "\n")
    print(f"Wrote {outname}")

    if args.profile:
        coords_files = sorted(glob.glob(f"{prefix}_*.coords"))
        if not coords_files:
            print(f"--profile requested but no {prefix}_*.coords files found "
                  "(the run needs -c/--coords-frequency set to produce them).")
        for path in coords_files:
            m = re.fullmatch(re.escape(prefix_basename) + r"_(.+)\.coords", os.path.basename(path))
            voltage = m.group(1) if m else "unknown"
            frames = read_coords_frames(path)
            centers, mean_hist = compute_profile(pore_obj, ptype, frames, args.profile_bins)
            profile_out = f"{prefix}_{voltage}_profile.dat"
            coord_label = "z[A]" if ptype == 'slit' else "r[A]"
            with open(profile_out, 'w') as f:
                f.write(f"# {'z-density' if ptype == 'slit' else 'radial-density'} "
                        f"profile from {path}\n")
                f.write(f"# Config: {os.path.abspath(args.config)}\n")
                f.write(f"# (1){coord_label}  (2)mean ion count per bin per frame\n")
                for c, h in zip(centers, mean_hist):
                    f.write(f"{c:.6e}  {h:.6e}\n")
            print(f"Wrote {profile_out}")


if __name__ == "__main__":
    main()
