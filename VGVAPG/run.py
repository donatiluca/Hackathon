from openmm import *
from openmm.app import *
from openmm.unit import *
from sys import stdout

import numpy as np
from time import gmtime, strftime

# directories
inp_dir  =  'input/'
out_dir  =  '/data/numerik/people/ldonati/OPENMM_ZIB/VGVAPG/first_test/'
out_dir  =  'C:/Users/donat/ZIB_DIR_TOO_BIG/OPENMM_ZIB/VGVAPG/first_test/'
out_dir  =  ''

# Integrator parameters
dt      = 0.002      # ps
Nsteps  = 2500000
Nout    = 50

print('*** Number of nanoseconds ***')
print(Nsteps * dt / 1000)

# System parameters
kB    = 0.008314  # kJ mol-1 K-1
T     = 300       # K
gamma = 1         # ps-1


# LOG-file
log = open(out_dir             + "log.txt", 'w')

log.write('Timestep: '         + str(dt)     + " ps\n")
log.write("nsteps: "           + str(Nsteps) + "\n" )
log.write("nstxout: "          + str(Nout)   + "\n")
log.write("Temperature: "      + str(T)      + " K\n")
log.write("Collision rate: "   + str(gamma)  + " ps-1\n")
log.write("Boltzmann const.: " + str(kB)     + " kJ mol-1 K-1 \n")
log.write("Simulation start: " + strftime("%Y-%m-%d %H:%M:%S", gmtime()) + "\n" )
log.close();

# input topology and gro file
prmtop = AmberPrmtopFile(inp_dir + 'VGVAPG_nowat.prmtop')
inpcrd = AmberInpcrdFile(inp_dir + 'VGVAPG_nowat.inpcrd')

# create system
system = prmtop.createSystem(implicitSolvent=OBC2, nonbondedMethod=CutoffNonPeriodic, nonbondedCutoff=1.0, constraints=None, rigidWater=True)


integrator = LangevinIntegrator(T*kelvin, gamma/picosecond, dt*picoseconds)

##########################################################################################################
###  S T A R T   S I M U L A T I O N  ####################################################################
##########################################################################################################

# set-up simulation
platform = Platform.getPlatformByName('Reference')
simulation = Simulation(prmtop.topology, system, integrator, platform)
simulation.context.setPositions(inpcrd.positions)

# minimization
print('\n\n*** Minimizing ...')
simulation.minimizeEnergy()
print('*** Minimization completed ***')

# equilibration
simulation.context.setVelocitiesToTemperature(T)
print('\n\n*** Equilibrating...')
simulation.step(10000)
print('*** Equilibration completed ***')




##########################################################################################################
###  S A V E S  ##########################################################################################
##########################################################################################################

# print on screen
simulation.reporters.append(StateDataReporter(stdout, 10000, speed = True, step=True, potentialEnergy=True, temperature=True))

# save trajectory
simulation.reporters.append(DCDReporter(out_dir + "trajectory.dcd", Nout))






# repeat procedure for nsteps
simulation.step(Nsteps)

# add total calculation time to LOG-file
log = open(out_dir + "log.txt", 'a')
log.write("Simulation end: " + strftime("%Y-%m-%d %H:%M:%S", gmtime()) )
log.close();

# end
print('\n\n****** SIMULATION COMPLETED *****************************\n\n')


