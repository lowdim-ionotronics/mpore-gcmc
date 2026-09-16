# Original code written by Taras Verkholyak,
# Institute for Condensed Matter Physics,
# National Academy of Sciences of Ukraine
#
# Adapted and modified by Vishnu Prasad Kurupath
# (https://vishnu-prasad-kurupath.github.io/)
#
# Voltage-sweep GCMC runner: runs a range of electrode potentials with a
# given step, carrying the final configuration from each voltage over as
# the initial configuration for the next. Automatically resumes from
# ./cont.restart if present. For a single-voltage run, see run_gcmc.py.

import sys
import os
import json
import numpy as np
import mpore_gcmc
from argparse import ArgumentParser

parser = ArgumentParser(fromfile_prefix_chars='@')
parser.add_argument("-o", "--output-prefix", help="output file prefix", metavar="FILE", type=str, dest="output_prefix")
parser.add_argument("-A", "--pore-width-accessible", help="accessible pore width in A -- diameter for --pore-type cyl, gap width for --pore-type slit -- i.e. the width available to ion centres", type=float, metavar="VAL", dest="pore_width_accessible")
parser.add_argument("--pore-width-nominal", help="nominal pore width in A, to the wall-atom/carbon centres (diameter for cyl, gap width for slit). Provide this and/or --pore-width-accessible; consistent with --wall-atom-radius or an error is raised", type=float, metavar="VAL", dest="pore_width_nominal")
parser.add_argument("-a", "--ion-radii", help="ion radii (in A)", type=float, metavar="A1, A2, ...", dest="ion_radii", nargs='+')
parser.add_argument("-w", "--mu-bulk", help="bulk chemical potential for each particle type (in eV)", type=float, nargs='+', metavar="MU1,MU2,...", dest="mu_bulk")
parser.add_argument("-u", "--voltage", help="voltage range start,stop,step (in V)", type=float, nargs='+', metavar="VAL", dest="voltage")
parser.add_argument("-n", "--production-steps", help="number of MC steps", metavar="VAL", type=int, dest="production_steps")
parser.add_argument("-t", "--thermalization-steps", help="number of MC steps for thermalization", metavar="VAL", type=int, dest="thermalization_steps")
parser.add_argument("-P", "--pore-type", help="pore type: cyl/slit", metavar="VAL", type=str, dest="pore_type")
parser.add_argument("-T", "--temperature", help="temp in K", metavar="VAL", type=float, dest="temperature")
parser.add_argument("-e", "--eshift", help="distance (in Angs) to shift electron center from pore atom center", metavar="VAL", type=float, dest="eshift")
parser.add_argument("-p", "--epsr", help="permittivity of medium", metavar="VAL", type=float, dest="epsr")
parser.add_argument("-L", "--pore-length", help="periodic simulation length in A -- cylinder axial length for --pore-type cyl, square lateral box side length for --pore-type slit", metavar="VAL", type=float, dest="pore_length")
parser.add_argument("-W", "--wall-atom-radius", help="radius of wall atoms in A (e.g. carbon radius); relates --pore-width-accessible and --pore-width-nominal", metavar="VAL", type=float, dest="wall_atom_radius")
parser.add_argument("-s", "--stat-frequency", help="how often to dump statistics", metavar="VAL", type=float, dest="stat_frequency")
parser.add_argument("-r", "--restart-frequency", help="how often to write restart", metavar="VAL", type=float, dest="restart_frequency")
parser.add_argument("-c", "--coords-frequency", help="how often to dump coordinates", metavar="VAL", type=float, dest="coords_frequency")
parser.add_argument("-pt", "--prob-trans", help="probability for MC translation moves", metavar="VAL", type=float, dest="prob_trans")
parser.add_argument("-pw", "--prob-widom", help="probability for MC widom moves", metavar="VAL", type=float, dest="prob_widom")
parser.add_argument("-C", "--config", help="JSON file of parameters (same names as the long options above); supplies defaults for any option not given on the command line", metavar="FILE", type=str, dest="config")
args = parser.parse_args()

if args.config:
    with open(args.config) as f:
        cfg = json.load(f)
    for key, value in cfg.items():
        if key.startswith('_'):
            continue
        if getattr(args, key, None) is None:
            setattr(args, key, value)

if args.output_prefix is not None:
    prefix = args.output_prefix
else:
    print("The prefix name (-o/--output-prefix) is not provided. Defaulting to \"Case\".")
    prefix = "Case"
if args.production_steps is not None:
    n_simul = args.production_steps
else:
    print("The number of simulation steps (-n/--production-steps) is not provided. Defaulting to 1000000.")
    n_simul = 1000000
if args.thermalization_steps is not None:
    n_therm = args.thermalization_steps
else:
    print("The number of thermalization steps (-t/--thermalization-steps) is not provided. Defaulting to 500000.")
    n_therm = 500000
if args.ion_radii:
    aion = np.array(args.ion_radii)
    if len(aion) < 2:
        print("Provide at least two ion radii A (-a/--ion-radia)")
        sys.exit()
    q = np.array([-1, 1])
    for i in np.arange(1, len(aion) - 1):
        q = np.append(q, 0)
    print("Charges: ", q)
else:
    print("Give the ion radia in A (-a/--ion-radii)")
    sys.exit()
if args.voltage:
    u = np.array(args.voltage)
else:
    print("Provide voltage in V (-u/--voltage)")
    sys.exit()
if args.mu_bulk:
    w = np.array(args.mu_bulk)
else:
    print("Provide bulk chemical potential in eV (-w/--mu-bulk)")
    sys.exit()
if args.pore_width_accessible is None and args.pore_width_nominal is None:
    print("Provide the pore width (-A/--pore-width-accessible and/or --pore-width-nominal)")
    sys.exit()
if args.pore_type is not None:
    ptype = args.pore_type
else:
    print("Provide the pore type (-P/--pore-type TYPE), where TYPE=slit/cyl")
    sys.exit()
if args.temperature is not None:
    temp = args.temperature
else:
    print("Provide the temperature (-T/--temperature)")
    sys.exit()
if args.eshift is not None:
    eshift = args.eshift
else:
    print("eshift (-e/--eshift) not provided. Defaulting to 0.0")
    eshift = 0.0
if args.epsr is not None:
    epsr = args.epsr
else:
    print("Permittivity (-p/--epsr) not provided. Defaulting to 2.5")
    epsr = 2.5
if args.pore_length is not None:
    Ltube = args.pore_length
else:
    print("Pore length (-L/--pore-length) not provided. Defaulting to 100.0")
    Ltube = 100.0
wall_atom_radius = args.wall_atom_radius if args.wall_atom_radius is not None else 0.0
try:
    pore_width_accessible, pore_width_nominal = mpore_gcmc.resolve_pore_width(
        args.pore_width_accessible, args.pore_width_nominal, wall_atom_radius)
except ValueError as e:
    print(e)
    sys.exit()
if args.stat_frequency is not None:
    stat_freq = args.stat_frequency
else:
    print("Statistics frequency (-s/--stat-frequency) not provided. Defaulting to 1000")
    stat_freq = 1000
if args.restart_frequency is not None:
    restart_freq = args.restart_frequency
else:
    print("Restart frequency (-r/--restart-frequency) not provided. Defaulting to 1000")
    restart_freq = 1000
if args.coords_frequency is not None:
    coord_dump = args.coords_frequency
else:
    print("Coord dump frequency (-c/--coords-frequency) not provided. Defaulting to None (no dump)")
    coord_dump = None
if args.prob_trans is not None:
    p_trans = args.prob_trans
else:
    print("Translation probability (-pt/--prob-trans) not provided. Defaulting to 0.5")
    p_trans = 0.5
if args.prob_widom is not None:
    p_widom = args.prob_widom
else:
    print("Widom probability (-pw/--prob-widom) not provided. Defaulting to 1.0")
    p_widom = 1.0

if u.size != 3:
    print("Voltage range incorrectly provided. Defaulting to [0.0, 0.8] with a step of 0.05")
    u = [0.0, 0.8, 0.05]
else:
    print("Voltage range", u)

state = mpore_gcmc.state(temp, epsr, aion, q,
                ptype, pore_width_accessible, Ltube, eshift, wall_atom_radius,
                n_therm, n_simul, p_trans, p_widom)

if os.path.isfile('./cont.restart'):
    state = mpore_gcmc.pickle_load('./cont.restart')

mc_exec = mpore_gcmc.mcfunctions(state)
del state

while mc_exec.state.c_voltage <= u[1]:

    if mc_exec.state.mcparams.complete == 0:

        ueV = mc_exec.state.c_voltage / mc_exec.state.factor_kBT_to_eV
        mu_comp = w / mc_exec.state.factor_kBT_to_eV + mc_exec.state.q_comp * ueV

        print('Voltage:', mc_exec.state.c_voltage)
        print('Total (mu+ev) for ions', mu_comp)

        print('Thermalization started from', mc_exec.state.mcparams.c_therm)
        mc_exec.thermalization(mu_comp, restart_freq, continuation=True)
        print('Thermalization ended')

        print('Production started from', mc_exec.state.mcparams.c_sim)
        mc_exec.mc_simulate(mu_comp, stat_freq, prefix, mc_exec.state.c_voltage, coord_dump, restart_freq, continuation=True)
        print('Production ended')

    mc_exec.state.c_voltage += u[2]
    mc_exec.state.mcparams.complete = 0
