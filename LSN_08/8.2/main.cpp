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
#include <iomanip>
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
   int n_grid=40;
   double mu_vect[n_grid];
   double sigma_vect[n_grid];
   double grid_energy[n_grid][n_grid][2];//Valor medio e incertezza nelle due componenti finali
   double center_mu=+0.8;
   double length_mu=+0.2;//Lunghezza totale dell'intervallo
   double center_sigma=+0.6;
   double length_sigma=+0.2;//Lunghezza totale dell'intervallo

   //Preparo la griglia
   for(int j=0;j<n_grid;j++){   
     mu_vect[j]=j*(length_mu/n_grid)+center_mu-length_mu/2;
     sigma_vect[j]=j*(length_sigma/n_grid)+center_sigma-length_sigma/2;
   }
   
   //Eseguo le simulazioni
   for(int j=0;j<n_grid;j++){
     
     if( (j*10)%n_grid==0){cout<<endl<<int((j*100)/n_grid)<<"% "<<endl;}
     
     for(int k=0;k<n_grid;k++){ 
       
       double mu=mu_vect[j];
       double sigma=sigma_vect[k];
       
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
	 
	 //if(i%(N/100)==0 and i<=(N/10)){cout<<endl<<"Accettazione: "<<conta_accetta/double(i+1)<<" a "<<i/(N/100)<<"% dei passi effettuati"<< endl;}
       }
       
       
       monteC stima(raw_H +dati_saltati,N-dati_saltati,100);//L'ultimo valore e' il numero dei blocchi
       
       stima.statistica();
       grid_energy[j][k][0]=stima.media();
       grid_energy[j][k][1]=stima.errore();
       //stima.stampa("Energy.dat");
       
     }
   }
   
   cout<<endl<<"100%"<<endl;
   
   //Ora trovo il minimo della griglia!
   
   int index_min[n_grid];//Vettore in cui raccogliere gli indici dei valori minimi a j fissato e k variabile

   //In questo ciclo cerco i minimi a j fissato e k variabile
   for(int j=0;j<n_grid;j++){
     index_min[j]=0;
     for(int k=0;k<n_grid;k++){
       if(grid_energy[j][k][0]<grid_energy[j][ index_min[j] ][0]){
	 index_min[j]=k;
       }
     }
   }

   //In questo invece faccio variare j tenendo i k racchiusi in index_min[]
   int index_min_2=0;
   for(int j=0;j<n_grid;j++){
     if(grid_energy[j][ index_min[j] ][0]<grid_energy[index_min_2][ index_min[index_min_2] ][0]){
       index_min_2=j;
     }
   }

   int index_sigma_opt=index_min[index_min_2];
   int index_mu_opt=index_min_2;

   cout<<endl<<"Minimum in the grid: mu = "<<mu_vect[index_mu_opt]<<", sigma = "<<sigma_vect[index_sigma_opt]<<", E = "<<grid_energy[index_mu_opt][index_sigma_opt][0]<<" ± "<<grid_energy[index_mu_opt][index_sigma_opt][1]<<endl;

   cout<<endl<<endl;

   const int wd=12;
   
   if(index_mu_opt!=0 and index_mu_opt != n_grid-1 and index_sigma_opt!=0 and index_sigma_opt != n_grid-1){
     cout<<"σ/μ"<<setw(wd)<<mu_vect[index_mu_opt-1]<<setw(wd)<<mu_vect[index_mu_opt]<<setw(wd)<<mu_vect[index_mu_opt+1]<<endl<<endl;
     cout<<sigma_vect[index_sigma_opt-1]<<setw(wd)<<grid_energy[index_mu_opt-1][index_sigma_opt-1][0]<<setw(wd)<<grid_energy[index_mu_opt][index_sigma_opt-1][0]<<setw(wd)<<grid_energy[index_mu_opt+1][index_sigma_opt-1][0]<<endl<<endl;
     cout<<sigma_vect[index_sigma_opt]<<setw(wd)<<grid_energy[index_mu_opt-1][index_sigma_opt][0]<<setw(wd)<<grid_energy[index_mu_opt][index_sigma_opt][0]<<setw(wd)<<grid_energy[index_mu_opt+1][index_sigma_opt][0]<<endl<<endl;
     cout<<sigma_vect[index_sigma_opt+1]<<setw(wd)<<grid_energy[index_mu_opt-1][index_sigma_opt+1][0]<<setw(wd)<<grid_energy[index_mu_opt][index_sigma_opt+1][0]<<setw(wd)<<grid_energy[index_mu_opt+1][index_sigma_opt+1][0]<<endl<<endl;
   }

   
   //Per finire rieseguo la simulazione sul dato più preciso e stampo le posizioni visitate.
   double mu=mu_vect[index_mu_opt];
   double sigma=sigma_vect[index_sigma_opt];
   
   unsigned int N=1000000; //Passi
   int dati_saltati=1000;
   N=N+dati_saltati;
   
   double x;
   double delta=3.;
   int conta_accetta=0;
   double raw_H[N];
   double* posizioni=new double[N];
   
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

     posizioni[i]=x;
     raw_H[i]=kinetic(x,mu,sigma)/psi(x,mu,sigma)+potential(x);
     
     //if(i%(N/100)==0 and i<=(N/10)){cout<<endl<<"Accettazione: "<<conta_accetta/double(i+1)<<" a "<<i/(N/100)<<"% dei passi effettuati"<< endl;}
   }
   
   
   monteC stima(raw_H +dati_saltati,N-dati_saltati,100);//L'ultimo valore e' il numero dei blocchi
   
   stima.statistica();
   stima.stampa("Energy.dat");

   cout<<endl<<"Minimum in the grid, second calculation: mu = "<<mu_vect[index_mu_opt]<<", sigma = "<<sigma_vect[index_sigma_opt]<<", E = "<<stima.media()<<" ± "<<stima.errore()<<endl;

   
   //Stampo su file le posizioni da dopo l'equilibrazione
   ofstream file_out;
   file_out.open("Psi^2.dat");
   for(unsigned int j=dati_saltati;j<N;j++){
     file_out<<posizioni[j]<<endl;
   }
   file_out.close();
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
