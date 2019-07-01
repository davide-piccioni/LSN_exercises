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

   
   unsigned int N=100000000; //Numero di lanci totali
   unsigned int k=100;//Numero di blocchi
   unsigned int L=N/k; //Elemeni per blocco

   //Definisco lunghezza dell'ago e spaziatura da usare nell'esperimento di Buffon
   const double width=0.5;//Lunghezza dell'ago
   const double d=1.;//Spaziatura

   //Ora eseguo N lanci e inserisco il risultato in un vettore (1 se incrocia una linea, 0 se non la incrocia).
   int * raw=new int[N];

   for(unsigned int j=0;j<N;j++){
     
     //Genero il centro di massa dell'ago tra 0 e d.
     double y_cm=rnd.Rannyu(0.5*d,1.5*d);

     //Potrei fare così, ma userei il pigreco per valutare pigreco
     //double theta=rnd.Rannyu(-M_PI/2,M_PI/2);

     //Quindi uso un metodo accept-reject
     double x=0.;
     double y=0.;
     do{
       x=rnd.Rannyu();
       y=rnd.Rannyu(-1,+1);
     }while(x*x+y*y>1 or (x==0 and y==0));//Genero due punti nella semi-circonferenza a destra di raggio unitario       
     
     double coseno=0.5*(x/sqrt(x*x+y*y));

     if(y_cm>d+width*coseno or y_cm<(d-width*coseno) ){
       raw[j]=0;
     }
     else{
       raw[j]=1;
     }

   }

   
   //Riempio i vettori che contengono le medie parziali e i quadrati delle stesse
   double * prog_1=new double[k];
   double * prog_2=new double[k];
   double * err_prog=new double[k];
   
   //Riempio prima il primo elemento
   int somma=0;
   for(unsigned int l=0;l<L;l++){
     somma+=raw[l];
     //cout<<somma<<endl;
   }
   prog_1[0]=(2*width*L)/(somma*d);//Inserisco la prima stima di pi greco
   prog_2[0]=(prog_1[0])*(prog_1[0]);//Ne inserisco il quadrato
   
   
     //Ora riempio tutti gli altri
     for(unsigned int i=1;i<k;i++){
       int sommadue=0;
       for(unsigned int j=0;j<(L);j++){
	 sommadue+=raw[i*L+j];
       }
       prog_1[i]=prog_1[i-1]+(2*width*L)/(sommadue*d);
       prog_2[i]=pow((2*width*L)/(sommadue*d),2) +prog_2[i-1];
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
