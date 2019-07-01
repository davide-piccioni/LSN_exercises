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
   
   unsigned int N=1000000; //Numero di numeri totali
   unsigned int k=100;//Numero di blocchi
   unsigned int L=N/k; //Elemeni per blocco
   
   //Riempio direttamente i vettori che contengono le medie parziali e i quadrati delle stesse
   double * prog_1=new double[k];
   double * prog_2=new double[k];
   double * err_prog=new double[k];
   
   //Riempio prima il primo elemento
   double somma=0;
     for(unsigned int l=0;l<L;l++){
       somma+=pow((rnd.Rannyu()-0.5),2);
     }
   prog_1[0]=somma/L;
   prog_2[0]=(prog_1[0])*(prog_1[0]);

   
   //Ora riempio tutti gli altri
   for(unsigned int i=1;i<k;i++){
     double sommadue=0;
       for(unsigned int j=0;j<(L);j++){
	 sommadue+=pow((rnd.Rannyu()-0.5),2);
       }
       prog_1[i]=(prog_1[i-1]+(sommadue/L));
       prog_2[i]=((sommadue/L)*(sommadue/L) +prog_2[i-1]);
   }
   
   //Ora trasformo le somme parziali in medie e stimo le incertezze
   for(unsigned int i=1;i<k;i++){
     prog_1[i]=prog_1[i]/(i+1);
     prog_2[i]=prog_2[i]/(i+1);
     err_prog[i]=sqrt((prog_2[i]-(prog_1[i])*(prog_1[i]))/(i));
   }
   err_prog[0]=0;

   //Stampo i dati
   for(unsigned int i=0;i<k;i++){
     file_out<<prog_1[i]<<' '<<err_prog[i]<<endl;
     //file_out<<err_prog[i]<<endl;
     //cout<<(err_prog[i])<<endl;
     //cout<<prog_1[i]<<endl;
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
