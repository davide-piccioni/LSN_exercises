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
#include <cstring>
#include <cmath>
#include "genetic.h"

using namespace std;

int main (){

  //Directory where to find input and send output
  //The input must always be of the form initial_cities.dat
  string locus="./Square_30/";

  //String operations so that changing the directory is implemented
  string _input=locus+"initial_cities.dat";
  char input[_input.size()+1];
  strcpy(input,&_input[0]);
  string _output=locus+"final_cities.dat";
  char output[_output.size()+1];
  strcpy(output,&_output[0]);
  string _L_best=locus+"L_best.dat";
  char L_best[_L_best.size()+1];
  strcpy(L_best,&_L_best[0]);
  string _L_average=locus+"L_average.dat";
  char L_average[_L_average.size()+1];
  strcpy(L_average,&_L_average[0]);
  
  gensalesman primo(input);
  int n_generations = 50;

  ofstream out_best;
  ofstream out_average;

  out_best.open(L_best);
  out_average.open(L_average);
  
  for(int generations=0;generations<n_generations;generations++){  
    primo.check_population();
    primo.sort_popul();
    
    cout<<endl<<"Generazione: "<<generations<<", L_best = "<<primo.best_L()<<endl;
    out_best<<primo.best_L()<<endl;
    out_average<<primo.average_L()<<endl;
    
    primo.marriage();
    primo.mutate();
    primo.check_population();
  }
  cout<<endl<<"Generazione: "<<n_generations<<", L_best = "<<primo.best_L()<<endl;
  primo.print_bestroute(output);
  out_best.close();
  out_average.close();
  
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
