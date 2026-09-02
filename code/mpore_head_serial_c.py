# Original code written by Taras Verkholyak, 
# Institute for Condensed Matter Physics, 
# National Academy of Sciences of Ukraine, 
# and Andrij Kuzmak, Department for Theoretical Physics, 
# Ivan Franko National University of Lviv
# 
# Adapted and modified by Vishnu Prasad Kurupath
# (https://vishnu-prasad-kurupath.github.io/)
 
import sys
import os
import numpy as np
import mpore_serial_c as mpore
from argparse import ArgumentParser

parser = ArgumentParser(fromfile_prefix_chars='@')
parser.add_argument("-o", "--output-prefix", help="output file prefix", metavar="FILE", type=str, dest="prefix")
parser.add_argument("-A", "--tube-radius", help="accessible tube radius (in A)", type=float, metavar="VAL", dest="atube")
parser.add_argument("-a", "--ion-radii", help="ion radia (in A)", type=float, metavar="A1, A2, ...", dest="aion", nargs='+')
parser.add_argument("-w", "--transfer-energy", help="transfer energy for each particle type (in eV)", type=float, nargs='+', metavar="W1,W2,...", dest="w")
parser.add_argument("-u", "--voltage", help="voltage range start,stop,step (in V)", type=float, nargs='+', metavar="VAL", dest="u")
parser.add_argument("-n", "--production-steps", help="number of MC steps", metavar="VAL", type=int, dest="n")
parser.add_argument("-t", "--thermalization-steps", help="number of MC steps for thermalization", metavar="VAL", type=int, dest="ntherm")
parser.add_argument("-P", "--pore-type", help="pore type: cyl/slit", metavar="VAL", type=str, dest="ptype")
parser.add_argument("-T", "--temperature", help="temp in K", metavar="VAL", type=float, dest="temp")
parser.add_argument("-e", "--eshift", help="distance (in Angs) to shift electron center from pore atom center", metavar="VAL", type=float, dest="eshift")
parser.add_argument("-p", "--epsr", help="permittivity of medium", metavar="VAL", type=float, dest="epsr")
parser.add_argument("-L", "--tube-length", help="length of tube (size size for slit)", metavar="VAL", type=float, dest="Ltube")
parser.add_argument("-W", "--wall-atom-radius", help="Radius of wall atoms to reduce from pore width", metavar="VAL", type=float, dest="wall_atom_radius")
parser.add_argument("-s", "--stat-frequency", help="How often to dump coordinates", metavar="VAL", type=float, dest="stat_freq")
parser.add_argument("-r", "--restart-frequency", help="How often to write restart", metavar="VAL", type=float, dest="restart_freq")
parser.add_argument("-c", "--coords-frequency", help="How often to dump coordinates", metavar="VAL", type=float, dest="coord_dump")
parser.add_argument("-pt", "--prob-trans", help="Probability for MC translation moves", metavar="VAL", type=float, dest="p_trans")
parser.add_argument("-pw", "--prob-widom", help="Probability for MC widom moves", metavar="VAL", type=float, dest="p_widom")
args = parser.parse_args()

if args.prefix:
    prefix = args.prefix
else:
    print ("The prefix name (-o/--output-prefix) is not provided. Defaulting to \"Case\".")
    prefix = "Case"
if args.n:
    n_simul = args.n
else:
    print ("The number of simulation steps (-n/--production-steps) is not provided. Defaulting to 1000000.")
    n_simul = 1000000 
if args.ntherm:
    n_therm = args.ntherm
else:
    print ("The number of thermalization steps (-t/--thermalization-steps) is not provided. Defaulting to 500000.")
    n_therm = 500000 
if args.aion:
    aion = np.array(args.aion)
    q = np.array([-1, 1])
    for i in np.arange(1,len(aion)-1):
        q=np.append(q,0);
    print ("Charges: ", q)
    if len(aion) < 2:
        print ("Provide at least two ion radii A (-a/--ion-radia)")
        sys.exit()
else:
    print ("Give the ion radia in A (-a/--ion-radii)")
    sys.exit()
if args.u:
    u = np.array(args.u)
else:
    print ("Provide voltage in V (-u/--voltage)")
    sys.exit()
if args.w:
    w = np.array(args.w)
else:
    print ("Provide transfer energy in KBT (-w/--transfer-energy)")
    sys.exit()
if args.atube:
    atube = args.atube
else:
    print ("Provide the tube radius in A (-A/--tube-radius)")
    sys.exit()
if args.ptype:
    ptype = args.ptype
else:
    print ("Provide the pore type (-P/--pore-type TYPE), where TYPE=slit/cyl")
    sys.exit()
if args.temp:
    temp = args.temp
else:
    print ("Provide the temperature (-T/--temperature)")
    sys.exit()
if args.eshift:
    eshift = args.eshift
else:
    print ("eshift (-e/--eshift) not provided. Defaulting to 0.0")
    eshift = 0.0
if args.epsr:
    epsr = args.epsr
else:
    print ("Permittivity (-p/--epsr) not provided. Defaulting to 2.5")
    epsr = 2.5
if args.Ltube:
    Ltube = args.Ltube
else:
    print ("Tube Length (-L/--tube-length) not provided. Defaulting to 100.0")
    Ltube = 100.0
if args.wall_atom_radius:
    wall_atom_radius = args.wall_atom_radius
else:
    print ("Wall atom radius (-W/--wall-atom-radius) not provided. Defaulting to 0.0")
    wall_atom_radius = 0.0
if args.stat_freq:
    stat_freq = args.stat_freq
else:
    print ("Stastics frequency (-s/--stat-frequency) not provided. Defaulting to 1000")
    stat_freq = 1000
if args.restart_freq:
    restart_freq = args.restart_freq
else:
    print ("Restart frequency (-r/--restart-frequency) not provided. Defaulting to 1000")
    restart_freq = 1000
if args.coord_dump:
    coord_dump = args.coord_dump
else:
    print ("Coord dump frequency (-c/--coords-frequency) not provided. Defaulting to None (no dump)")
    coord_dump = None
if args.p_trans:
    p_trans = args.p_trans
else:
    print ("Translation probability (-pt/--prob-trans) not provided. Defaulting to 0.5")
    p_trans = 0.5
if args.p_widom:
    p_widom = args.p_widom
else:
    print ("Widom probability (-pw/--prob-widom) not provided. Defaulting to 1.0")
    p_widom = 1.0

if u.size != 3:
    print("Voltage range incorrectly provided. Defaulting to [0.0, 0.8] with a step of 0.05")
    u = [0.0, 0.8, 0.05]
else:
    print("Voltage range", u)

state = mpore.state(temp, epsr, aion, q, 
                ptype, 2*atube, Ltube, eshift, wall_atom_radius, 
                n_therm, n_simul, p_trans, p_widom)

if os.path.isfile('./cont.restart'):
    state = mpore.pickle_load()

mc_exec = mpore.mcfunctions(state)
del state

while mc_exec.state.c_voltage <= u[1]:
    
    if mc_exec.state.mcparams.complete == 0:

        ueV = mc_exec.state.c_voltage/mc_exec.state.factor_kBT_to_eV
        mu_comp = w/mc_exec.state.factor_kBT_to_eV+mc_exec.state.q_comp*ueV

        print('Voltage:', mc_exec.state.c_voltage)
        print('Total (mu+ev) for ions', mu_comp)

        print('Thermalization started from', mc_exec.state.mcparams.c_therm)
        mc_exec.thermalization(mu_comp, restart_freq)
        print('Thermalization ended')

        print('Production started from', mc_exec.state.mcparams.c_sim)
        mc_exec.mc_simulate(mu_comp, stat_freq, prefix, mc_exec.state.c_voltage, coord_dump, restart_freq)
        print('Production ended')

    mc_exec.state.c_voltage += u[2]
    mc_exec.state.mcparams.complete = 0
