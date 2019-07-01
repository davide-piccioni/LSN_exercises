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
#include <string>
#include <cmath>
#include "monteC.h"


using namespace std;

double psi(double x, double mu, double sigma);
double psi_2(double x, double mu, double sigma);
double kinetic(double x, double mu, double sigma);
double potential(double x);

int main (int argc, char *argv[]){
  
  //Inizializzo generatore
  ////////////////////////
  Random rnd;
  int seed[4];
  int p1, p2;
  ifstream Primes("Primes");
  if (Primes.is_open()){
    Primes >> p1 >> p2 ;
  } else cerr << "PROBLEM: Unable to open Primes" << endl;
  Primes.close();

  ifstream input("seed.in");
  string property;
  if (input.is_open()){
    while ( !input.eof() ){
      input >> property;
      if( property == "RANDOMSEED" ){
	input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
	rnd.SetRandom(seed,p1,p2);
      }
    }
    input.close();
  } else cerr << "PROBLEM: Unable to open seed.in" << endl;
  //////////////////////////////////////
  
  //Parametri
  double mu=0.8;
  double sigma=0.6;
  
  unsigned int N=1000000; //Passi
  int dati_saltati=1000;
  N=N+dati_saltati;
  
  double x;
  double delta=3.;
  int conta_accetta=0;
  double raw_H[N];
  
  //Posizione iniziale
  x=0.;
  
  for(unsigned int i=0; i<N; i++){//Passi
    
    //if(i%(N/100)==0){cout<<i<<" su "<<N<<endl;}
    
    //Genero x uniformemente distribuito con passo delta
    double x_new;
    
    //x_new=rnd.Gauss(x,delta);
    x_new=x+delta*rnd.Rannyu(-1,+1);
    
    //Uso Metropolis
    double alpha=psi_2(x_new,mu,sigma)/psi_2(x,mu,sigma);
    //if(alpha>1){alpha=1;}//Serve a fare il minimo, ma e' inutile.
    if(rnd.Rannyu()<=alpha){
      x=x_new;
      conta_accetta++;
    }
    
    raw_H[i]=kinetic(x,mu,sigma)/psi(x,mu,sigma)+potential(x);
    
    if(i%(N/100)==0 and i<=(N/10)){cout<<endl<<"Accettazione: "<<conta_accetta/double(i+1)<<" a "<<i/(N/100)<<"% dei passi effettuati"<< endl;}
  }
  
  
  monteC stima(raw_H+dati_saltati,N-dati_saltati,100);//L'ultimo valore e' il numero dei blocchi
  
  stima.statistica();
  
  stima.stampa("Energy.dat");
  
  
  return 0;
}

double psi(double x, double mu, double sigma){
  return exp( -pow((x-mu),2)/(2*pow(sigma,2)) )+exp( -pow((x+mu),2)/(2*pow(sigma,2)) );
}

double psi_2(double x, double mu, double sigma){
  return exp( -pow((x-mu),2)/pow(sigma,2) )+exp( -pow((x+mu),2)/pow(sigma,2) )+2*exp( -x*x-mu*mu/pow(sigma,2) );
}

double kinetic(double x, double mu, double sigma){
  return (1/(2*sigma*sigma))*(  exp( -pow((x-mu),2)/(2*pow(sigma,2)) )*( 1-pow((x-mu),2)/pow(sigma,2) ) + exp( -pow((x+mu),2)/(2*pow(sigma,2)) )*( 1-pow((x+mu),2)/pow(sigma,2) )   );
}

double potential(double x){
  return pow(x,4)-2.5*pow(x,2);
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
