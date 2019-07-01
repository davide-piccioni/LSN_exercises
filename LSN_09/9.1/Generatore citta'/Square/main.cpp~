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
#include "random.h"


using namespace std;
 
int main (int argc, char *argv[]){

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

   /////////////////////////////////////////////////
   //Da qui inizia il programma
   ofstream file_out;

   file_out.open("square_cities.dat");
   
   unsigned int N=30; //Numero di citta'
   
   file_out<<N<<endl;
   
   double x[N];
   double y[N];
   
   for(unsigned int l=0;l<N;l++){
       x[l]=rnd.Rannyu(-1,+1);
       y[l]=rnd.Rannyu(-1,+1);
     }     
     
   //Stampo i dati
   for(unsigned int i=0;i<N;i++){
     file_out<<x[i]<<' '<<y[i]<<endl;
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
