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
   ///////////////////////////////////////
   const unsigned int N=100000; //Ripetizioni esperimento

   const double T=1; //Tempo
   const double mu=0.1; //Risk-free interest rate
   const double sigma=0.25; //Volatilita'
   const double S_0=100; //Prezzo al tempo zero
   const unsigned int M=100; //Passi in cui discretizzare l'intervallo

   double* raw_call=new double[N];//Vettore in cui inserire il risultato della singola simulazione per opzione call
   double* raw_put=new double[N];//Vettore in cui inserire il risultato della singola simulazione per opzione put

   monteC stima_call(N,100);//Il secondo valore e' il numero dei blocchi
   monteC stima_put(N,100);//Il secondo valore e' il numero dei blocchi

   //Genero i risultati degli esperimenti
     for(unsigned int l=0;l<N;l++){
       //Genero il dato
       double S_T=S_0;

       for(unsigned int j=0;j<M;j++){
	 S_T=S_T*exp( (mu-0.5*sigma*sigma)*(T/M) + sigma*rnd.Gauss(0,1)*sqrt(T/M));
       }

       //Opzioni call e put
       if( (S_T-S_0)>0){
	 raw_call[l]=exp(-mu*T)*(S_T-S_0);
	 raw_put[l]=0.;
       }
       else{
	 raw_put[l]=exp(-mu*T)*(S_0-S_T);
	 raw_call[l]=0.;
       }
     }

     //Ora uso gli oggetti per stimare le incertezze con il metodo dei blocchi
     stima_call.set_raw(raw_call);
     stima_put.set_raw(raw_put);
     stima_call.statistica();
     stima_put.statistica();

     cout<<"Call: "<<stima_call.media()<<" "<<stima_call.errore()<<endl;
     cout<<"Put: "<<stima_put.media()<<" "<<stima_put.errore()<<endl;

     //Stampa i valori intermedi dei blocchi
     stima_call.stampa("risultati_call.dat");
     stima_put.stampa("risultati_put.dat");

   return 0;
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
