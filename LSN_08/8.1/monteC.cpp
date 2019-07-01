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
#include <cmath>
#include <cstdlib>
#include "monteC.h"


using namespace std;

monteC :: monteC(){
  
  ////////////////////////////////////
  //Inizializzo generatore di numeri casuali
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
  ////////////////////////////////////

  N=100000;
  k=100;
  L=N/k;
  raw=new double[N];
  prog_1=new double[k];
  prog_2=new double[k];
  err_prog=new double[k];
}

monteC :: monteC(unsigned int N_tuo, unsigned int k_tuo){
  
  ////////////////////////////////////
  //Inizializzo generatore di numeri casuali
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
  ////////////////////////////////////

  N=N_tuo;
  k=k_tuo;
  L=N/k;
  raw=new double[N];
  prog_1=new double[k];
  prog_2=new double[k];
  err_prog=new double[k];
}

monteC :: monteC(double * my_raw, unsigned int N_tuo, unsigned int k_tuo){
  
  ////////////////////////////////////
  //Inizializzo generatore di numeri casuali
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
  ////////////////////////////////////

  N=N_tuo;
  k=k_tuo;
  L=N/k;
  raw=new double[N];
  prog_1=new double[k];
  prog_2=new double[k];
  err_prog=new double[k];

  for(unsigned int j=0;j<N;j++){
    raw[j]=my_raw[j];
}

}

monteC :: monteC(char * input_file, unsigned int N_tuo, unsigned int k_tuo){
  
  ////////////////////////////////////
  //Inizializzo generatore di numeri casuali
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
  ////////////////////////////////////

  N=N_tuo;
  k=k_tuo;
  L=N/k;
  raw=new double[N];
  prog_1=new double[k];
  prog_2=new double[k];
  err_prog=new double[k];

  ifstream file_in;
  file_in.open(input_file);
  if(file_in.fail()){
    cout<<"Problema apertura file "<< input_file <<endl;
    return;
  }

  for(unsigned int j=0;j<N;j++){
    file_in>>raw[j];
}
  file_in.close();
}


monteC :: ~monteC(){
  delete [] raw;
  delete [] prog_1;
  delete [] prog_2;
  delete [] err_prog;
}

double** monteC :: statistica(){

  //Sommo i primi L elementi
  double somma_prima=0;
  for(unsigned int j=0;j<L;j++){
    somma_prima+=raw[j];
  }
  
  //Riempio prima il primo elemento
  prog_1[0]=somma_prima/L;
  prog_2[0]=(prog_1[0])*(prog_1[0]);//Ne inserisco il quadrato
  
  
  //Ora riempio tutti gli altri
  for(unsigned int i=1;i<k;i++){
    double somma=0;
    for(unsigned int j=0;j<L;j++){
      somma+=raw[L*i+j];
    }
    prog_1[i]=prog_1[i-1]+somma/L;
    prog_2[i]=(somma/L)*(somma/L) +prog_2[i-1];
  }
  
  //Ora trasformo le somme parziali in medie e stimo le incertezze
  for(unsigned int i=1;i<k;i++){
    prog_1[i]=prog_1[i]/(i+1);
    prog_2[i]=prog_2[i]/(i+1);
    err_prog[i]=sqrt((prog_2[i]-(prog_1[i])*(prog_1[i]))/(i));
  }
  err_prog[0]=0;

  double** risultati=new double*[2];
  risultati[0]=new double[k];
  risultati[1]=new double[k];
  
    for(unsigned int j=0;j<k;j++){
      risultati[0][j]=prog_1[j];
    }
  for(unsigned int j=0;j<k;j++){
    risultati[1][j]=err_prog[j];
  }

  return risultati;  
}

void monteC :: stampa_schermo()const {
  for(unsigned int i=0;i<k;i++){
    cout<<prog_1[i]<<' '<<err_prog[i]<<endl;
  }
  return;
}

void monteC :: stampa(const char* output)const {
  ofstream file_out;
  
  file_out.open(output);
  
  for(unsigned int i=0;i<k;i++){
    file_out<<prog_1[i]<<' '<<err_prog[i]<<endl;
  }
  
  file_out.close();
  return;
}

void monteC ::set_raw (double * my_raw){
  for(unsigned int j=0;j<N;j++){
    raw[j]=my_raw[j];
}
  return;
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
