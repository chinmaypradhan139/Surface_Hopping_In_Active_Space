
Description:
This repository contains two folders.
1. Surface_Hopping_in_Active_Space
   The folder Surface_Hopping_in_Active_Space contains all the files necessary to reproduce the population data of Surface_Hopping_In_Active_Space, as well as the Marcus theory.
   The input file for both Surface_Hopping_In_Active_Space and Marcus theory is AFSSH.inp. The input file is self-explanatory.
   To run the surface hopping simulation, choose appropriate parallelization in the first line of the input file, then type ./run.
   When the jobs finish, open AFSSH.inp select the first line as 3 and run ./run to average the populations.
   This will generate a file pop.out containing the population of spin-up and spin-down states with time.
   Matrix_evolve.py is a short script to get the population of four different Marcus states. On running this file, depending on the input, you can get the results of populations of different states as in the paper.
	 mod_afssh.f90 contains all the subroutines to run surface hopping.
	 
2. IESH_2025
  	fort.23 is the input file where all quantities have to be input in atomic units.
    Running this script would output the population vs time file in each folder as fort . 104.
    To run this, use ./myscript.sh
    To average the population, simply run the average.py
   
   



Credits:
I, Chinmay Pradhan, have developed the code for Surface hopping in active space using the AFSSH code, which was created by Amber Jain and Aarti Sindhu in Prof. Amber Jain's group at IIT Bombay.
The link to that AFSSH code is https://github.com/amber-jain-group-iitb/AFSSH-7site-anharmonic-bath.
