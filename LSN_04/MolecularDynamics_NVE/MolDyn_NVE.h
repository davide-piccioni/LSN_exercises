/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
#ifndef __MolDyn_NVE__
#define __MolDyn_NVE__

const int nbins=100;

class MolDyn{

private:

//parameters, observables
  double stima_pot, stima_kin, stima_etot, stima_temp, stima_press;

/*
//bins for g(r)
int bins[nbins];
double bin_size;
*/

// averages
double acc,att;

//configuration
double* x, *y,*z,*xold,*yold,*zold;
double* vx,*vy,*vz;

// thermodynamical state
int npart;
double energy,temp,vol,rho,box,rcut;

// simulation
int nstep, iprint;
double delta;

public:

 double get_temp()const{return temp;}
 int get_npart()const{return npart;}
 double get_rho()const{return rho;}
 double get_rcut()const{return rcut;}
 double get_box()const{return box;}
 double get_delta()const{return delta;}
 int get_nstep()const{return nstep;}
 int get_iprint()const{return iprint;}
 
 MolDyn(const char* input, const char* config);//Constructor from file for the first time -> input: parameters; config: positions
//ACHTUNG: Random velocities are generated according to Maxwell-Blotzmann at inserted temperature
 
 MolDyn(const char* input,const char* config1, const char* config2);//Constructor from file -> input: parameters; config1, config2: positions
 //ACHTUNG: Gets last two positions of the system and restarts from there rescaling velocities.
 
 void Move(void);//Algorithm for moving
 void Measure(void);//Algorithm for measuring
 double Force(int ip, int idir);//Function to calculate force on particle ip over direction idir
 double Pbc(double);
 
 void ConfFinal(const char* output);//Print final positions
 void ConfFinal_and_previous(const char* output1, const char* output2);//Print final positions and previous one);//Print final positions and previous one
 void ConfXYZ(int);
};


#endif

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
