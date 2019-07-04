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
   ofstream file_out;
   file_out.open("risultati.dat");

   unsigned int N=1000000; //Ripetizioni esperimento
   unsigned int k=100; //Passi per ogni catena

   double * x=new double[N];//Coordinata x, ripetizione i-esima
   double * y=new double[N];//Idem
   double * z=new double[N];//Idem

   double* raw=new double[N];//Distanza per tutte le catene a passo fissato

   monteC stima(N,100);//Il secondo valore e' il numero dei blocchi

   //Inizializzo tutti
   for(unsigned int l=0;l<N;l++){
     x[l]=0.;
     y[l]=0.;
     z[l]=0.;
     raw[l]=0.;
   }

   file_out<<0.<<' '<<0.<<' '<<0.<<endl;

   for(unsigned int i=0; i<k; i++){//Passi

     cout<<"Passo "<<i<<endl;

     for(unsigned int l=0;l<N;l++){//Catena

       double theta=acos(1-2*rnd.Rannyu());
       double phi=rnd.Rannyu(0,2*M_PI);

       z[l]+=cos(theta);
       x[l]+=sin(theta)*cos(phi);
       x[l]+=sin(theta)*sin(phi);

     //Ora faccio la distanza per ogni catena al passo k
       raw[l]=(x[l]*x[l]+y[l]*y[l]+z[l]*z[l]); 
       //cout<<"Catena "<<l<<": r**2= "<<raw[l]<<endl;          
     }
     
     //cout<<"Debug"<<endl;

     //Ora creo un oggetto che per ogni passo stimi le incertezze con il metodo dei blocchi
     stima.set_raw(raw);
     stima.statistica();
     double media_quad=stima.media();
     double err_quad=stima.errore();
     //file_out<<sqrt(media_quad)<<' '<<1/(2*sqrt(media_quad))*err_quad<<endl;
     //cout<<sqrt(media_quad)<<' '<<1/(2*sqrt(media_quad))*err_quad<<endl;
     file_out<<sqrt(media_quad)<<' '<<1/(2*sqrt(media_quad))*err_quad<<' '<<1/(2*media_quad)*err_quad<<endl;
     //stima.stampa_schermo();
     //Propagazione degli errori, ale'!
   }

   file_out.close();
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
