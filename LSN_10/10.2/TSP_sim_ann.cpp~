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
#include <algorithm> //Random_shuffle, sort
#include "TSP_sim_ann.h"


anneal_tsp :: anneal_tsp(const char* file_input){
  
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

  ifstream file_in;
  file_in.open(file_input);
  
  if(file_in.fail() ){
  cout<<"Errore in lettura file"<<endl;
  }
  
  file_in>>N;
  
  our_path.alleles.resize(N);

  cities_x.resize(N);
  cities_y.resize(N);
  
  //Leggo le posizioni delle citta'
  for(int i=0;i<N;i++){
    file_in>>cities_x[i];
    file_in>>cities_y[i];
  }
  
  //Genero il percorso casualmente
  for(int i=0;i<N;i++){
    our_path.alleles[i]=i+1;
  }
  random_shuffle ( our_path.alleles.begin(), our_path.alleles.end() );
  
  /*
  cout<<endl;
  for(int j=0;j<N;j++){
    cout<<our_path.alleles[j]<<' ';
  }
  cout<<endl;
  */

  calculate_L();
  
}


anneal_tsp :: anneal_tsp(const char* file_input, int n_riga){
  
  ////////////////////////////////////
  //Inizializzo generatore di numeri casuali
  int seed[4];
  int p1, p2;
  ifstream Primes("Primes");
  if (Primes.is_open()){
    for(int i=0;i<n_riga;i++){
      Primes >> p1 >> p2 ;
    }
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

  ifstream file_in;
  file_in.open(file_input);
  
  if(file_in.fail() ){
  cout<<"Errore in lettura file"<<endl;
  }
  
  file_in>>N;
  
  our_path.alleles.resize(N);

  cities_x.resize(N);
  cities_y.resize(N);
  
  //Leggo le posizioni delle citta'
  for(int i=0;i<N;i++){
    file_in>>cities_x[i];
    file_in>>cities_y[i];
  }
  
  //Genero il percorso casualmente
  for(int i=0;i<N;i++){
    our_path.alleles[i]=i+1;
  }

  for(int l=0;l<2*N;l++){
    int indice_rand_1=rnd.Rannyu()*N;
    int indice_rand_2=rnd.Rannyu()*N;
    int appoggio=our_path.alleles[indice_rand_1];
    our_path.alleles[indice_rand_1]=our_path.alleles[indice_rand_2];
    our_path.alleles[indice_rand_2]=appoggio;
  }
  
  /*
  cout<<endl;
  for(int j=0;j<N;j++){
    cout<<our_path.alleles[j]<<' ';
  }
  cout<<endl;
  */

  calculate_L();

}

void anneal_tsp :: calculate_L(){

 //Ora calcolo L(x)
  our_path.L=0.;
  //Calcola L tra tutte le città contigue nel vettore
  for(int j=0;j<N-1;j++){
    our_path.L += sqrt(  pow(cities_x[ our_path.alleles[j+1]-1 ]-cities_x[ our_path.alleles[j]-1 ],2) + pow(cities_y[ our_path.alleles[j+1]-1 ]-cities_y[ our_path.alleles[j]-1 ],2)  );
  }
  //Calcola L anche tra la prima e l'ultima citta'
  our_path.L += sqrt(  pow(cities_x[ our_path.alleles[N-1]-1 ]-cities_x[ our_path.alleles[0]-1 ],2) + pow(cities_y[ our_path.alleles[N-1]-1 ]-cities_y[ our_path.alleles[0]-1 ],2)  );
  //cout<<endl<<our_path.L<<endl;

  return;
}

double anneal_tsp :: get_L(){
  calculate_L();
  return our_path.L;
}

bool anneal_tsp :: check_path(){

  bool check_param=0;

  vector <int> prova(N);
  
  prova=our_path.alleles;
  sort( prova.begin(), prova.end() );
  for(int j=0;j<N;j++){
    if(prova[j]!=j+1 and check_param==0){
      check_param=1;
    }
  }

  if(check_param==1){
    cout<<endl<<"Problema con il percorso: valori non sensati!!!"<<endl;
  }
  
  return check_param;
}

bool anneal_tsp :: exchange(double beta){

  calculate_L();

  pathway proposed_path=our_path;

  int mutated_1=rnd.Rannyu()*N;
  int mutated_2=0;
  do{mutated_2=rnd.Rannyu()*N;}while(mutated_2==mutated_1);    
  int appo= proposed_path.alleles[mutated_1];
  proposed_path.alleles[mutated_1]=proposed_path.alleles[mutated_2];
  proposed_path.alleles[mutated_2]=appo;

  //Ora calcolo L(x) per il nuovo percorso
  proposed_path.L=0.;
  //Calcola L tra tutte le città contigue nel vettore
  for(int j=0;j<N-1;j++){
    proposed_path.L += sqrt(  pow(cities_x[ proposed_path.alleles[j+1]-1 ]-cities_x[ proposed_path.alleles[j]-1 ],2) + pow(cities_y[ proposed_path.alleles[j+1]-1 ]-cities_y[ proposed_path.alleles[j]-1 ],2)  );
  }
  //Calcola L anche tra la prima e l'ultima citta'
  proposed_path.L += sqrt(  pow(cities_x[ proposed_path.alleles[N-1]-1 ]-cities_x[ proposed_path.alleles[0]-1 ],2) + pow(cities_y[ proposed_path.alleles[N-1]-1 ]-cities_y[ proposed_path.alleles[0]-1 ],2)  );

  //Uso Metropolis
  double alpha=exp(-beta*proposed_path.L)/exp(-beta*our_path.L);
  if(rnd.Rannyu()<=alpha){
    our_path=proposed_path;
    //Mossa accettata
    return 1;
  }

  //In questo caso mossa non accettata
  return 0;
}

bool anneal_tsp :: local_shift(double beta){
 
  calculate_L();

  pathway proposed_path=our_path;

  int first_index=int( rnd.Rannyu()*(N-2) );
  int last_index=int( rnd.Rannyu(1,(N-first_index)/2) );
  rotate( proposed_path.alleles.begin()+first_index,proposed_path.alleles.begin()+first_index+last_index,proposed_path.alleles.begin() + first_index + 2*last_index );
  
  //Ora calcolo L(x) per il nuovo percorso
  proposed_path.L=0.;
  //Calcola L tra tutte le città contigue nel vettore
  for(int j=0;j<N-1;j++){
    proposed_path.L += sqrt(  pow(cities_x[ proposed_path.alleles[j+1]-1 ]-cities_x[ proposed_path.alleles[j]-1 ],2) + pow(cities_y[ proposed_path.alleles[j+1]-1 ]-cities_y[ proposed_path.alleles[j]-1 ],2)  );
  }
  //Calcola L anche tra la prima e l'ultima citta'
  proposed_path.L += sqrt(  pow(cities_x[ proposed_path.alleles[N-1]-1 ]-cities_x[ proposed_path.alleles[0]-1 ],2) + pow(cities_y[ proposed_path.alleles[N-1]-1 ]-cities_y[ proposed_path.alleles[0]-1 ],2)  );

  //Uso Metropolis
  double alpha=exp(-beta*proposed_path.L)/exp(-beta*our_path.L);
  if(rnd.Rannyu()<=alpha){
    our_path=proposed_path;
    //Mossa accettata
    return 1;
  }

  //In questo caso mossa non accettata
  return 0;

}

bool anneal_tsp :: local_permut(double beta){

  calculate_L();

  pathway proposed_path=our_path;

  int first_index=int( rnd.Rannyu(0,N-3) );
  int second_index=int( rnd.Rannyu(first_index+1,N-2)  );
  int third_index=int( rnd.Rannyu(second_index+1,N-1) );
  int last_index=int( rnd.Rannyu(third_index+1,N) );
  
  vector <int> appoggio((1+second_index-first_index)+(1+last_index-third_index));
  int conta=0;
  for(int k=first_index;k<=second_index;k++){
    appoggio[conta]=proposed_path.alleles[k];
    conta++;
  }
  for(int k=third_index;k<=last_index;k++){
    appoggio[conta]=proposed_path.alleles[k];
    conta++;
  }
  
  random_shuffle( appoggio.begin(),appoggio.end() );
  
  conta=0;
  for(int k=first_index;k<=second_index;k++){
    proposed_path.alleles[k]=appoggio[conta];
    conta++;
  }
  for(int k=third_index;k<=last_index;k++){
    proposed_path.alleles[k]=appoggio[conta];
    conta++;
  }
  
  //Ora calcolo L(x) per il nuovo percorso
  proposed_path.L=0.;
  //Calcola L tra tutte le città contigue nel vettore
  for(int j=0;j<N-1;j++){
    proposed_path.L += sqrt(  pow(cities_x[ proposed_path.alleles[j+1]-1 ]-cities_x[ proposed_path.alleles[j]-1 ],2) + pow(cities_y[ proposed_path.alleles[j+1]-1 ]-cities_y[ proposed_path.alleles[j]-1 ],2)  );
  }
  //Calcola L anche tra la prima e l'ultima citta'
  proposed_path.L += sqrt(  pow(cities_x[ proposed_path.alleles[N-1]-1 ]-cities_x[ proposed_path.alleles[0]-1 ],2) + pow(cities_y[ proposed_path.alleles[N-1]-1 ]-cities_y[ proposed_path.alleles[0]-1 ],2)  );

  //Uso Metropolis
  double alpha=exp(-beta*(proposed_path.L-our_path.L) );
  if(rnd.Rannyu()<=alpha){
    our_path=proposed_path;
    //Mossa accettata
    return 1;
  }

  //In questo caso mossa non accettata
  return 0;

}


vector <double> anneal_tsp :: move(double beta, int n_steps){

  vector <double> acceptance_rate(3);
  acceptance_rate[0]=0.;
  acceptance_rate[1]=0.;
  acceptance_rate[2]=0.;

  for(int j=0;j<n_steps;j++){
    double die=rnd.Rannyu();
    if(die<1./3){
      acceptance_rate[0]+=exchange(beta);
    }
    else{
      if(die<2./3){
	acceptance_rate[1]+=local_shift(beta);
      }
      else{
	acceptance_rate[2]+=local_permut(beta);
      }
    } 
  }

  acceptance_rate[0]/=n_steps;
  acceptance_rate[1]/=n_steps;
  acceptance_rate[2]/=n_steps;

  return acceptance_rate;
}


void anneal_tsp :: print_path_numbers(){
  
  cout<<endl;
  for(int j=0;j<N;j++){
    cout<<our_path.alleles[j]<<' ';
  }
  cout<<endl;

  return;
}

void  anneal_tsp :: print_route(const char* outfile){
  
  ofstream file_out;
  file_out.open(outfile);
  for(int l=0;l<N;l++){
    file_out<<cities_x[our_path.alleles[l]-1]<<" "<<cities_y[our_path.alleles[l]-1]<<endl;
  }
}

anneal_tsp :: ~anneal_tsp(){

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
