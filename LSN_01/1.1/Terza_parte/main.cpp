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

   file_out.open("risultati.dat");
   
   unsigned int M=100;//Numero di bin
   unsigned int L=10000; //Elementi per ogni chi
   unsigned int K=100;//Numero test chi-quadro da fare

   double * chi=new double[K];

   for(unsigned int j=0;j<K;j++){//Singolo Chi-quadro
     
     //Questo è il vettore dei bin, a cui aggiungo un 1 per ogni numero nel bin corrispondente.
     int * bin=new int[M];
     for(unsigned int m=0;m<M;m++){//Inizializza tutti i valori
       bin[m]=0;
     }
     
     for(unsigned int l=0;l<L;l++){//L lanci
       double lancio=rnd.Rannyu();
       bin[int(lancio*M)]+=1;//Riempio i bin in questo modo
     }

     //Ora calcolo il chi-quadro
     for(unsigned int m=0;m<M;m++){//Inizializza tutti i valori
       chi[j]+=pow( (bin[m]-double(L/M)), 2)/(L/M);
     }
     
   }
   

   //Stampo i dati
   for(unsigned int i=0;i<K;i++){
     file_out<<chi[i]<<endl;
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
