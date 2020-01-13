using namespace std;

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <iostream>
#include <fstream>
#include <istream>
#include <list>
#include <ostream>
#include <math.h>
#include <vector>
#include <random>			//for random generator
#include <functional>		//for random generator as function
#include <Eigen/Dense>	//library "Eigen" necessary to use the matrices/vectors for the Langevin propagation

int main(int argc, char* argv[]){
    
for (int i = 1; i < argc; ++i) {
  if(std::string(argv[i]) == "-h"){
    cout << " " <<endl;
    cout << "USAGE: ./LE_1dim_reflect [OPTIONS]" << endl;
    cout << " " <<endl;
    cout << "INPUT: files containing negative logarithm of the histogram representing the free energy (in kJ/mol), friction constant (in kJ ps /(mol ns**2 nm**2))," << endl;
    cout << "start point (in nm), mass (in kg/mol), temperature (in K). Free energy and friction units are compatible with input from dcTMD correction scripts" << endl; 
    cout << "available at www.moldyn.uni-freiburg.de/software/software.html." << endl; 
    cout << "OUTPUT: trajectory generated according to the Markovian Langevin equation (in ns, nm)" <<endl;
    cout << "RESOLUTION: integration timestep (option -t) times output frequency (option -s) is the time difference between the points." <<endl;
    cout << " " <<endl;
    cout << "The program propagates the Markovian Langevin equation in one dimension based on given fields starting at a certain point (file at option -start)." << endl;
    cout << " " <<endl;
    cout << "The Markovian Langevin equation is based on the free energy (drift), Stokes friction and white noise. Based on F(x)=-kT*ln(P(x)), the file given to the program" <<endl;
    cout << "(option -free) needs to be (x,kJ/mol). The program needs in addition the number of free energy and friction points in the respective files" <<endl;
    cout << "(specified by options -n and -ngamma). x needs to be given in nm." <<endl;
    cout << "The derivative of the free energy is approximated by the difference between the local free energy and the free energy at the neighbouring grid points." <<endl;
    cout << "The differences in both directions (x_i > x_{i,actual} and x_i < x_{i,actual}) are averaged. The friction gamma (units of kJ ps /(mol ns**2 nm**2))" <<endl;
    cout << "is given by the file at option -gamma similar to the free energy, the noise is scaled according to the Fluctuation-Dissipation theorem." <<endl;
    cout << "Together with the temperature (option -T given in K), the mass (option -mass given in kg/mol) and the integration timestep (option -t given in ns)" <<endl;
    cout << "this is everything which we need. Based on normal distributed white noise (seed given by option -I), the velocity is propagated by according to the integration" <<endl;
    cout << "scheme presented by Bussi and Parrinello in Bussi G. and Parrinello M., 'Accurate sampling using Langevin dynamics', Phys. Rev. E, 75(5), 056707, (2007)." <<endl;
    cout << " " <<endl;
    cout << "The constructed trajectory is reflected at the borders of the free energy file. If the trajectory jumps over one of these borders by the" <<endl;
    cout << "distance a, it is set back to x_{shift}(t)=x(t)-a. In addition, the velocity is multiplied by (-1) to mimick an elastic collision." <<endl;
    cout << " " <<endl;
    cout << " " <<endl;
    cout << "ATTENTION: The friction (given by -gamma) needs as additional information the number of bins (-ngamma) so that it can deviate from the free energy in this regard." <<endl;
    cout << "ATTENTION: If the numbers given with -n and -ngamma do not fit to the files, the code will not notice this problem. So, pay attention." <<endl;
    cout << " " <<endl;
    cout << " " <<endl;
    cout << "OPTIONS: after each option there needs to be a whitespace, e.g., '-T 320', to not disturb the program" <<endl;
    cout << " " << " " << "-h show these lines" << endl;
    cout << " " << " " << "-start name file with starting point [default startpoint in nm]" << endl;
    cout << " " << " " << "-free name file with free energy (x,G(x)) [default free_energy, kJ/mol]" << endl;
    cout << " " << " " << "-gamma name file with friction (x,Gamma(x)) [default gamma, units kg/(mol*ns), equals the kJ ps /(mol nm^2) from Gromacs]" << endl;
    cout << " " << " " << "-mass name file with mass [default mass, units kg/mol]" << endl;
    cout << " " << " " << "-o name output trajectory [default lang, units ns,nm]" << endl;
    cout << " " << " " << "-t integration timestep in ns [default 1 fs, i.e., 0.000001]" << endl;
    cout << " " << " " << "-T temperature [default 300 K]" << endl;
    cout << " " << " " << "-I seed random number generator [default: 0]" << endl;
    cout << " " << " " << "-L length of output trajectory [default: 100000 points]" << endl;
    cout << " " << " " << "-s write out every sth point [default: 1, i.e., every timestep]" << endl;
    cout << " " << " " << "-n number of 'grid' points of the free energy [default: 200]" << endl;
    cout << " " << " " << "-ngamma number of 'grid' points of gamma [default: 200]" << endl;\
    cout << " " <<endl;
    cout << " " <<endl;
    exit(0);
  }
}

// SCRIPT ASSUMES LENGTH SCALES AND TIMES AS NM AND NS, IN THIS WAY IT RESCALES EVERYTHING

//Count Variables
int a,c,f;

//Input Variables
int length,freeenergybins,gammabins,everynth;
int seed=0;
double deltat=0.0;
double temperature=0.0;
    
std::string namefreeenergy="free_energy";
std::string namegamma="gamma";
std::string namemass="mass";
std::string nameoutput="lang";
std::string namestart="startpoint";

length=100000;
seed=0;
deltat=0.000001;
temperature=300;
freeenergybins=200;
gammabins=200;
everynth=1;

for(a=1;a<argc;a++){
    if(std::string(argv[a])=="-free"){
        if(a+1< argc){ 							// Make sure we aren't at the end of argv!
		namefreeenergy=argv[a+1]; 						// Increment 'i' so we don't get the argument as the next argv[i].
        }
    }
    else if(std::string(argv[a])=="-o"){
        if(a+1< argc){ 
		nameoutput=argv[a+1];
        }
    }
    else if(std::string(argv[a])=="-start"){
        if(a+1< argc){ 
		namestart=argv[a+1];
        }
    }
    else if(std::string(argv[a])=="-gamma"){
        if(a+1< argc){ 
		namegamma=argv[a+1];
        }
    }
    else if(std::string(argv[a])=="-mass"){
        if(a+1< argc){ 
		namemass=argv[a+1];
        }
    }
    else if(std::string(argv[a])=="-max"){
        if(a+1<argc){
            seed=(int)(atof(argv[a+1]));
        }
    }
    else if(std::string(argv[a])=="-L"){
        if(a+1<argc){
		length=(int)(atof(argv[a+1]));
        }
    }
    else if(std::string(argv[a])=="-s"){
        if(a+1<argc){
		everynth=(int)(atof(argv[a+1]));
        }
    }
    else if(std::string(argv[a])=="-I"){
        if(a+1<argc){
		seed=(int)(atof(argv[a+1]));
        }
    }
    else if(std::string(argv[a])=="-t"){
        if(a+1<argc){
		deltat=atof(argv[a+1]);
        }
    }
    else if(std::string(argv[a])=="-T"){
        if(a+1<argc){
		temperature=atof(argv[a+1]);
        }
    }
    else if(std::string(argv[a])=="-n"){
        if(a+1<argc){
		freeenergybins=(int)(atof(argv[a+1]));
        }
    }
    else if(std::string(argv[a])=="-ngamma"){
        if(a+1<argc){
		gammabins=(int)(atof(argv[a+1]));
        }
    }
}

//lengthpropagate=length*everynth;
int numberfreeenergypoints=freeenergybins;
int numbergammapoints=gammabins;
double kT=8.3144*temperature;      // kT now in J/mol = kg*nm^2/(ns^2*mol)
    
//Input fields

double transfer=0.0;

double startpoint=0.0;
double maxcoordinates=0.0;
double mincoordinates=0.0;

double deltax=0.0;
double maxcoordinatesgamma=0.0;
double mincoordinatesgamma=0.0;

double deltaxgamma=0.0;
double massesinvers=0.0;
    
Eigen::VectorXd gammavary(numbergammapoints);
Eigen::VectorXd freeenergy(numberfreeenergypoints);
Eigen::VectorXd drift(numberfreeenergypoints);
Eigen::VectorXd C1(numberfreeenergypoints);
Eigen::VectorXd C2(numberfreeenergypoints);

ifstream datastart;                                                             //open input file
datastart.open(namestart,ios::in);
ifstream datamass;                                                             //open input file
datamass.open(namemass,ios::in);

datastart >> transfer;
startpoint=transfer;
datamass >> transfer;
massesinvers=1/transfer;

datastart.close();
datamass.close();

ifstream datafreeenergy;                                                             //open input file
datafreeenergy.open(namefreeenergy,ios::in);

for(a=0;a<numberfreeenergypoints;a++){
    datafreeenergy >> transfer;
    if(a==0){
        mincoordinates=transfer;
    }
    else if(a==numberfreeenergypoints-1){
        maxcoordinates=transfer;
    }
    datafreeenergy >> transfer;
    freeenergy(a)=transfer;
}

datafreeenergy.close();

ifstream datagamma;                                                             //open input file
datagamma.open(namegamma,ios::in);

for(a=0;a<numbergammapoints;a++){
    datagamma >> transfer;
    if(a==0){
        mincoordinatesgamma=transfer;
    }
    else if(a==numbergammapoints-1){
        maxcoordinatesgamma=transfer;
    }
    datagamma >> transfer;
    gammavary(a)=transfer;
}

datagamma.close();

if(maxcoordinates<=mincoordinates){
    cout << "Problems with the free energy format. Please check.";
    exit(0);
}
if(maxcoordinatesgamma<=mincoordinatesgamma){
    cout << "Problems with the gamma format. Please check.";
    exit(0);
}

mincoordinates=mincoordinates-0.000001;
maxcoordinates=maxcoordinates+0.000001;
deltax=(maxcoordinates-mincoordinates)/freeenergybins;
mincoordinatesgamma=mincoordinatesgamma-0.000001;
maxcoordinatesgamma=maxcoordinatesgamma+0.000001;
deltaxgamma=(maxcoordinatesgamma-mincoordinatesgamma)/gammabins;

cout << "Starting point, Free energy, gamma, tau read in." << endl;
cout << "Starting point:" << endl;
cout << startpoint << " ";
cout << endl;
cout << "Free energy range:" << endl;
cout << "Min:" << endl;
cout << mincoordinates << " ";
cout << endl;
cout << "Max:" << endl;
cout << maxcoordinates << " ";
cout << endl;
cout << "Delta x:" << endl;
cout << deltax << " ";
cout << endl;
cout << "Gamma range:" << endl;
cout << "Min:" << endl;
cout << mincoordinatesgamma << " ";
cout << endl;
cout << "Max:" << endl;
cout << maxcoordinatesgamma << " ";
cout << endl;
cout << "Delta x:" << endl;
cout << deltaxgamma << " ";
cout << endl;
cout << "Massinvers:" << endl;
cout << massesinvers << " ";
cout << endl;


//precalculation of drift term
cout << "calculating drift" << endl;

double actualfreeenergyplusdeltax=0.0,actualfreenergyminusdeltax=0.0;

for(a=1;a<(numberfreeenergypoints-1);a++){
    actualfreenergyminusdeltax=freeenergy(a-1);
    actualfreeenergyplusdeltax=freeenergy(a+1);
    drift(a)=-((actualfreeenergyplusdeltax-actualfreenergyminusdeltax)/deltax)*0.5*1000.;         // PMF gradient calculation; factor of 1000 for tansfer into J/mol
}
drift(0)=drift(1);
drift(numberfreeenergypoints-1)=drift(numberfreeenergypoints-2);

cout << "calculating propagation coefficients" << endl;
for(a=0;a<numberfreeenergypoints;a++){
    C1(a) = exp(-gammavary(a)*massesinvers*deltat*0.5); //in native Gromacs units
    C2(a) = sqrt((1-C1(a)*C1(a))*kT/massesinvers);
}


//random number generator

if(seed<0){                                                                     //mt19337 works just with unsigned ints, I don't know exactly the influences of negative seeds
    seed=-seed;
}

std::function<float()> randomnumber;							                 //random generator as function
randomnumber=std::bind(std::normal_distribution<double>(0.0, 1.0)			//uniform distributions of real numbers between zero and one
                 , std::mt19937(seed));

cout << "Start propagation" << endl;

ofstream output;                                                    //open output file
output.open(nameoutput,ios::out | ios::trunc);

output << "# t    x1" << endl;

int actualpointbinnumbers=0;
double pastpoint=0.0;
double actualpoint=0.0;
double pastvelocity=0.0;
double actualvelocity=0.0;
double pastdrift=0.0;
double actualdrift=0.0;
double pastwhitenoise=0.0;
double actualwhitenoise=0.0;
double pastC1=0.0;
double actualC1=0.0;
double pastC2=0.0;
double actualC2=0.0;

actualpoint=startpoint;

c=0;
int alreadymentioned=0;
double counter=0.0;
double timepoint=0.0;

for(a=0;a<length;a++){
    for(f=0;f<everynth;f++){
        if(c%100000==0){
            if(alreadymentioned==0){
                cout << "LE2 point" << " " << c << endl;
                alreadymentioned++;
            }
        }
        if(f==0){
            output << timepoint << " " << actualpoint;
            output << endl;
            c++;
            alreadymentioned=0;
        }
        actualdrift=drift(actualpointbinnumbers); 
        pastC1=C1(actualpointbinnumbers); // Gamma 
        pastC2=C2(actualpointbinnumbers); // Gamma 
        pastwhitenoise=randomnumber();
        actualwhitenoise=randomnumber();

                            // Bussi-integrator from Bussi G. and Parrinello M., 'Accurate sampling using Langevin dynamics', Phys. Rev. E, 75(5), 056707, (2007)
        pastpoint=actualpoint;
        pastvelocity=actualvelocity;
        pastdrift=actualdrift;

        actualpoint= pastpoint + (pastC1*pastvelocity + pastC2*massesinvers*pastwhitenoise)*deltat + massesinvers*pastdrift*deltat*deltat*0.5;

                            // reflection at box borders the box is defined by the upper and lower borders of the free energy file given to the program
        if(actualpoint>maxcoordinates){
            actualpoint=maxcoordinates-(actualpoint-maxcoordinates);
            pastvelocity=(-1)*pastvelocity;
        }
        if(actualpoint<mincoordinates){
            actualpoint=mincoordinates+(mincoordinates-actualpoint);
            pastvelocity=(-1)*pastvelocity;
        }

        actualpointbinnumbers=(int)((actualpoint-mincoordinates)/(maxcoordinates-mincoordinates)*freeenergybins);
        actualC1=C1(actualpointbinnumbers); // Gamma 
        actualC2=C2(actualpointbinnumbers); // Gamma
        actualdrift=drift(actualpointbinnumbers);
        actualvelocity=  actualC1*((pastC1*pastvelocity + pastC2*massesinvers*pastwhitenoise) + massesinvers*(pastdrift+actualdrift)*deltat*0.5) + actualC2*massesinvers*actualwhitenoise;

        counter = counter + 1.0;
        timepoint=counter*deltat;
    }
}

output.close();
return EXIT_SUCCESS;
}
