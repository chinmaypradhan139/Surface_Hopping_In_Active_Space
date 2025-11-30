
Repository Description

This repository contains two modules:

1. Surface_Hopping_in_Active_Space

This directory includes all files required to reproduce the population dynamics reported for Surface Hopping in Active Space, as well as the Marcus theory results.

Input file: AFSSH.inp (self-explanatory).

Running the simulations:
• Set the desired parallelization in the first line of AFSSH.inp, then execute:

./run


• After all jobs complete, set the first line of AFSSH.inp to 3 and run ./run again to compute the averaged populations.
• This will generate pop.out, which contains the time-resolved population of the spin-up and spin-down states.

Marcus theory states:
The script Matrix_evolve.py computes the populations of the four different Marcus states. Running the script with the desired input configuration will reproduce the population data shown in the manuscript.

Codebase:
The file mod_afssh.f90 contains the subroutines required to execute the surface hopping simulations.

2. IESH_2025

Input file: fort.23 (all parameters must be specified in atomic units).

Running the simulations:

./myscript.sh


This will produce the population-versus-time file fort.104 in each simulation directory.

Population averaging:

python3 average.py

The output.txt will contain the averaged population of the molecular electronic state.

Credits

The implementation of Surface Hopping in Active Space was developed by Chinmay Pradhan, building on the AFSSH code created by Amber Jain and Aarti Sindhu in Prof. Amber Jain’s group, IIT Bombay.
The original AFSSH repository is available at:
https://github.com/amber-jain-group-iitb/AFSSH-7site-anharmonic-bath
