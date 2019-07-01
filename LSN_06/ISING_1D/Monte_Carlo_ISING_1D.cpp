/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <iomanip>
#include "Monte_Carlo_ISING_1D.h"

using namespace std;

int main()
{
  //////////////////////////////////////////
  const int wd=12;
  for(double l=0.5;l<=2.01;l=l+0.01){
    ofstream file_out;
    file_out.open("input.dat");
    file_out<<l<<endl;//Temperatura
    file_out<<50<<endl;//nspin
    file_out<<1.0<<endl;//J
    file_out<<0.<<endl;//h
    file_out<<1<<endl;//1 per Metropolis, 0 per Gibbs
    file_out<<20<<endl;//nblk
    file_out<<10000<<endl;//nstep (Steps for each block)
    file_out.close();
  //////////////////////////////////////////
  
  Input("config.final"); //Inizialization: insert "0" for generating random spins, "config.something" to get it from that specific file
  
  for(int iblk=1; iblk <= nblk; ++iblk) //Simulation
  {
    cout<<endl<<"Temperatura = "<<l<<endl;
    Reset(iblk);   //Reset block averages
    for(int istep=1; istep <= nstep; ++istep)
    {
      Move(metro);
      Measure();
      Accumulate(); //Update block averages
    }
    Averages(iblk);   //Print results for current block
  }
  ConfFinal(); //Write final configuration

  ofstream out_t;
  out_t.open("T-U.dat",ios::app);
  out_t<<l<<setw(wd)<< glob_av[iu]/(double)nblk << setw(wd) << err_u << endl;
  out_t.close();
  out_t.open("T-C.dat",ios::app);
  out_t<<l<<setw(wd)<< glob_av[ic]/(double)nblk << setw(wd) << err_c << endl;
  out_t.close();
  out_t.open("T-M.dat",ios::app);
  out_t<<l<<setw(wd)<< glob_av[im]/(double)nblk << setw(wd) << err_m << endl;
  out_t.close();
  if(h==0.){
  out_t.open("T-X.dat",ios::app);
  out_t<<l<<setw(wd)<< glob_av[ix]/(double)nblk << setw(wd) << err_x << endl;
  out_t.close();
  }
  
  //////
  }
  //////
  return 0;
}


void Input(char* previous)//Se previous inserisci "0", genera nuovi spin random, altrimenti prendi dal file previous
{
  ifstream ReadInput;
  
  cout<<endl;
  cout << "Classic 1D Ising model             " << endl;
  cout << "Monte Carlo simulation             " << endl << endl;
  cout << "Nearest neighbour interaction      " << endl << endl;
  cout << "Boltzmann weight exp(- beta * H ), beta = 1/T " << endl << endl;
  cout << "The program uses k_B=1 and mu_B=1 units " << endl;

//Read seed for random numbers
   int p1, p2;
   ifstream Primes("Primes");
   Primes >> p1 >> p2 ;
   Primes.close();

   ifstream input("seed.in");
   input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
   rnd.SetRandom(seed,p1,p2);
   input.close();
  
//Read input informations
  ReadInput.open("input.dat");

  ReadInput >> temp;
  beta = 1.0/temp;
  cout << "Temperature = " << temp << endl;

  ReadInput >> nspin;
  cout << "Number of spins = " << nspin << endl;

  ReadInput >> J;
  cout << "Exchange interaction = " << J << endl;

  ReadInput >> h;
  cout << "External field = " << h << endl << endl;
    
  ReadInput >> metro; // if=1 Metropolis else Gibbs

  ReadInput >> nblk;

  ReadInput >> nstep;//Steps for each block

  if(metro==1) cout << "The program performs Metropolis moves" << endl;
  else cout << "The program perform Gibbs moves" << endl;
  cout << "Number of blocks = " << nblk << endl;
  cout << "Number of steps in one block = " << nstep << endl << endl;
  ReadInput.close();


//Prepare arrays for measurements
  iu = 0; //Energy
  iu2= 1; //Energy squared
  ic = 2; //Heat capacity
  im = 3; //Magnetization
  ix = 4; //Magnetic susceptibility
 
  n_props = 5; //Number of observables

  if(previous=="0"){
//initial configuration
    cout<<endl<<"Starting a new random configuration."<<endl;
    for (int i=0; i<nspin; ++i)
      {
	if(rnd.Rannyu() >= 0.5) s[i] = 1;
	else s[i] = -1;
      }
  }
  else{
    cout<<endl<<"Restarting from "<<previous<<" configuration."<<endl;
    ifstream file_in;
    file_in.open(previous);
    for (int i=0; i<nspin; ++i)
      {
	file_in>>s[i];
      }
    file_in.close();
  }
  
//Evaluate energy etc. of the initial configuration
  Measure();

//Print initial values for the potential energy and virial
  cout << "Initial energy = " << walker[iu]/(double)nspin << endl;
}


void Move(int metro)
{
  int o;

  for(int i=0; i<nspin; ++i)//Monte Carlo step!!!
    {
      //Select randomly a particle (for C++ syntax, 0 <= o <= nspin-1)
      o = (int)(rnd.Rannyu()*nspin);
      
      if(metro==1) //Metropolis
	{
	  double s_prop=0.;
	  if(s[o]==1.){s_prop=-1;}else{s_prop=+1.;}
	  double alpha=exp(-beta*(Boltzmann(s_prop,o)-Boltzmann(s[o],o)));
	  //if(alpha>1){alpha=1;}
	  if(rnd.Rannyu()<=alpha){
	    s[o]=s_prop;
	    accepted++;
	  }
	  attempted++;
	}
      else //Gibbs sampling
	{
	  //Alpha è la probabilità che o sia +1.
	  double alpha=1/(1+exp(+2*beta*(Boltzmann(+1,o)) ) );
	  if(rnd.Rannyu()<=alpha){
	    s[o]=+1;
	  }
	  else{
	    s[o]=-1;
	  }
	}
    }
}

double Boltzmann(int sm, int ip)
{
  double ene = -J * sm * ( s[Pbc(ip-1)] + s[Pbc(ip+1)] ) - h * sm;
  return ene;
}

void Measure()
{
  //int bin;
  double u = 0.0,  m = 0.0;

//cycle over spins
  for (int i=0; i<nspin; ++i)
  {
     u += -J * s[i] * s[Pbc(i+1)] - 0.5 * h * (s[i] + s[Pbc(i+1)]);
     m+=s[i];
  }
  walker[iu] = u;
  walker[iu2]=u*u;
  walker[im]=m;
  walker[ix]=m*m;
// INCLUDE YOUR CODE HERE
}


void Reset(int iblk) //Reset block averages
{
   
   if(iblk == 1)
   {
       for(int i=0; i<n_props; ++i)
       {
           glob_av[i] = 0;
           glob_av2[i] = 0;
       }
   }

   for(int i=0; i<n_props; ++i)
   {
     blk_av[i] = 0;
   }
   blk_norm = 0;
   attempted = 0;
   accepted = 0;
}


void Accumulate(void) //Update block averages
{

  for(int i=0; i<n_props; ++i)
    {
      blk_av[i] = blk_av[i] + walker[i];
    }
  blk_norm = blk_norm + 1.0;
}


void Averages(int iblk) //Print results for current block
{
    
   ofstream Ene;
   const int wd=12;
    
    cout << "Block number " << iblk << endl;
    if(metro==1){cout << "Acceptance rate " << accepted/attempted << endl << endl;}
    
    Ene.open("output.ene.0",ios::app);
    stima_u = blk_av[iu]/blk_norm/(double)nspin; //Energy
    glob_av[iu]  += stima_u;
    glob_av2[iu] += stima_u*stima_u;
    err_u=Error(glob_av[iu],glob_av2[iu],iblk);
    Ene << setw(wd) << iblk <<  setw(wd) << stima_u << setw(wd) << glob_av[iu]/(double)iblk << setw(wd) << err_u << endl;
    Ene.close();

    Ene.open("output.heat.0",ios::app);
    stima_u2 = blk_av[iu2]/blk_norm/((double)nspin); //Energy squared
    glob_av[iu2]  += stima_u2;
    glob_av2[iu2] += stima_u2*stima_u2;
    stima_c=beta*beta*((glob_av[iu2]/(double)iblk)-((double)nspin *glob_av[iu]/(double)iblk)*(glob_av[iu]/(double)iblk));
    glob_av[ic]  += stima_c;
    glob_av2[ic] += stima_c*stima_c;
    err_c=Error(glob_av[ic],glob_av2[ic],iblk);
    Ene << setw(wd) << iblk <<  setw(wd) << stima_c << setw(wd) << glob_av[ic]/(double)iblk  << setw(wd) << err_c << endl;
    Ene.close();

    Ene.open("output.magn.0",ios::app);
    stima_m = blk_av[im]/blk_norm/((double)nspin); //Magnetization
    glob_av[im]  += stima_m;
    glob_av2[im] += stima_m*stima_m;
    err_m=Error(glob_av[im],glob_av2[im],iblk);
    Ene << setw(wd) << iblk <<  setw(wd) << stima_m << setw(wd) << glob_av[im]/(double)iblk << setw(wd) << err_m << endl;
    Ene.close();

    if(h==0.){
    Ene.open("output.susc.0",ios::app);
    stima_x = beta*blk_av[ix]/blk_norm/((double)nspin); //Magnetic susceptibility
    glob_av[ix]  += stima_x;
    glob_av2[ix] += stima_x*stima_x;
    err_x=Error(glob_av[ix],glob_av2[ix],iblk);
    Ene << setw(wd) << iblk <<  setw(wd) << stima_x << setw(wd) << glob_av[ix]/(double)iblk << setw(wd) << err_x << endl;
    Ene.close();
}

// INCLUDE YOUR CODE HERE

    cout << "----------------------------" << endl << endl;
}


void ConfFinal(void)
{
  ofstream WriteConf;

  cout << "Print final configuration to file config.final " << endl << endl;
  WriteConf.open("config.final");
  for (int i=0; i<nspin; ++i)
  {
    WriteConf << s[i] << endl;
  }
  WriteConf.close();

  rnd.SaveSeed();
}

int Pbc(int i)  //Algorithm for periodic boundary conditions
{
    if(i >= nspin) i = i - nspin;
    else if(i < 0) i = i + nspin;
    return i;
}

double Error(double sum, double sum2, int iblk)
{
    return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)iblk);
}

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
