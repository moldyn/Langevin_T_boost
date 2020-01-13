# Description of relevant files

1) error_estimation_T_boosting.ipynb: This jupyter notebook allows for the estimation
of the error at 300 K when using the concept of T-boosting.
2) LE_1dim_reflect.cpp: implemetation of the one-dimensional Markovian Langevin 
equation with the possibility to have a elastic reflection at a lower and a higher border
3) CMakeLists.txt: used to compile the .cpp code


# Quick Start

The jupyter notebook can be used out of the box.

To compile the C++ code, make sure that the Eigen-C++-library is installed on your system!

Then, perform the following steps:\
0) Download the C++ file together with CMakeLists.txt to some folder Langevin
1) Run "cmake -H. -Bbuild" in the shell
2) Run "cmake --build build -- -j3" in the shell

After this you get the executable LE_1dim_reflect inside of Langevin/build.

# Use of the program LE_1dim_reflect

By running "./LE_1dim_reflect -h" you get some explanation how to use the program inside your shell.

For the impatient users, the following text shows what "./LE_1dim_reflect -h" gives:

USAGE: ./LE_1dim_reflect [OPTIONS]

INPUT: files containing negative logarithm of the histogram representing the free energy (in kJ/mol), friction constant (in kJ ps /(mol ns^2 nm^2)), start point (in nm), mass (in kg/mol), temperature (in K). Free energy and friction units are compatible with input from dcTMD correction scripts available at www.moldyn.uni-freiburg.de/software/software.html.

OUTPUT: trajectory generated according to the Markovian Langevin equation (in ns, nm).

RESOLUTION: integration timestep (option -t) times output frequency (option -s) is the time difference between the points.

The program propagates the Markovian Langevin equation in one dimension based on given fields starting at a certain point (file at option -start).
     
The Markovian Langevin equation is based on the free energy (drift), Stokes friction and white noise. Based on F(x)=-kT*ln(P(x)), the file given to the program (option -free) needs to be (x,kJ/mol). The program needs in addition the number of free energy and friction points in the respective files (specified by options -n and -ngamma). x needs to be given in nm.
The derivative of the free energy is approximated by the difference between the local free energy and the free energy at the neighbouring grid points.
The differences in both directions (x_i > x_{i,actual} and x_i < x_{i,actual}) are averaged. The friction gamma (units of kJ ps /(mol ns^2 nm^2)) is given by the file at option -gamma similar to the free energy, the noise is scaled according to the Fluctuation-Dissipation theorem.
Together with the temperature (option -T given in K), the mass (option -mass given in kg/mol) and the integration timestep (option -t given in ns) this is everything which we need. Based on normal distributed white noise (seed given by option -I), the velocity is propagated by according to the integration scheme presented by Bussi and Parrinello in Bussi G. and Parrinello M., 'Accurate sampling using Langevin dynamics', Phys. Rev. E, 75(5), 056707, (2007).
     
The constructed trajectory is reflected at the borders of the free energy file. If the trajectory jumps over one of these borders by the
distance a, it is set back to x_{shift}(t)=x(t)-a. In addition, the velocity is multiplied by (-1) to mimick an elastic collision.
     
     
ATTENTION: The friction (given by -gamma) needs as additional information the number of bins (-ngamma) so that it can deviate from the free energy in this regard.

ATTENTION: If the numbers given with -n and -ngamma do not fit to the files, the code will not notice this problem. So, pay attention.
     
     
OPTIONS: after each option there needs to be a whitespace, e.g., '-T 320', to not disturb the program

-h show these lines\
-start name file with starting point [default startpoint in nm]\
-free name file with free energy (x,G(x)) [default free_energy, kJ/mol]\
-gamma name file with friction (x,Gamma(x)) [default gamma, units kg/(mol*ns), equals the kJ ps /(mol nm^2) from Gromacs]\
-mass name file with mass [default mass, units kg/mol]\
-o name output trajectory [default lang, units ns,nm]\
-t integration timestep in ns [default 1 fs, i.e., 0.000001]\
-T temperature [default 300 K]\
-I seed random number generator [default: 0]\
-L length of output trajectory [default: 100000 points]\
-s write out every sth point [default: 1, i.e., every timestep]\
-n number of 'grid' points of the free energy [default: 200]\
-ngamma number of 'grid' points of gamma [default: 200]
