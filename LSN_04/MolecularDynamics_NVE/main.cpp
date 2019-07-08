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
#include "MolDyn_NVE.h"
#include "monteC.h"

using namespace std;

int main(){

  //Directory where to find input and send output
  string locus="./Solid/";

  //String operations so that changing the directory is implemented
  string _input=locus+"input.dat";
  char input[_input.size()+1];
  strcpy(input,&_input[0]);
  string _config_f0=locus+"config.final.0";
  string _config_f1=locus+"config.final.1";
  char config_f0[_config_f0.size()+1];
  char config_f1[_config_f1.size()+1];
  strcpy(config_f0,&_config_f0[0]);
  strcpy(config_f1,&_config_f1[0]);
  
  /*
  //For starting only with particles positions and create random velocities
  string _config_i=locus+"config.0";
  char config_i[_config_i.size()+1];
  strcpy(config_i,&_config_i[0]);
  
  MolDyn primo(input,config_i);//Inizialization 
  int nconf = 1;
  for(int istep=1; istep <= primo.get_nstep(); ++istep){
     primo.Move();           //Move particles with Verlet algorithm
     if(istep%primo.get_iprint() == 0) cout << "Number of time-steps: " << istep << endl;
     if(istep%10 == 0){
        primo.Measure();     //Properties measurement
        //primo.ConfXYZ(nconf);//Write actual configuration in XYZ format //Commented to avoid "filesystem full"! 
        nconf += 1;
     }
  }
  */
  
  //For starting from two previous particle positions
  MolDyn primo(input,config_f0,config_f1);//Inizialization 
  int nconf = 1;
  for(int istep=1; istep <= primo.get_nstep(); ++istep){
     primo.Move();           //Move particles with Verlet algorithm
     if(istep%primo.get_iprint() == 0) cout << "Number of time-steps: " << istep << endl;
     if(istep%10 == 0){
        primo.Measure();     //Properties measurement
        //primo.ConfXYZ(nconf);//Write actual configuration in XYZ format //Commented to avoid "filesystem full"! 
        nconf += 1;
     }
  }
  

  //For calculating mean values
  const int blocks=40;
  
  monteC ekin("output_ekin.dat",primo.get_nstep()/10,blocks);//Input, elementi, numero di blocchi
  monteC epot("output_epot.dat",primo.get_nstep()/10,blocks);//Input, elementi, numero di blocchi
  monteC etot("output_etot.dat",primo.get_nstep()/10,blocks);//Input, elementi, numero di blocchi
  monteC temp("output_temp.dat",primo.get_nstep()/10,blocks);//Input, elementi, numero di blocchi
  monteC press("output_press.dat",primo.get_nstep()/10,blocks);//Input, elementi, numero di blocchi

  ekin.statistica();
  epot.statistica();
  etot.statistica();
  temp.statistica();
  press.statistica();
  
  cout<<endl;
  cout<<"Kinetic Energy: "<<ekin.media()<<" "<<ekin.errore()<<endl;
  cout<<"Potential Energy: "<<epot.media()<<" "<<epot.errore()<<endl;
  cout<<"Total Energy: "<<etot.media()<<" "<<etot.errore()<<endl;
  cout<<"Temperature: "<<temp.media()<<" "<<temp.errore()<<endl;
  cout<<"Pressure: "<<press.media()<<" "<<press.errore()<<endl;

  double a;
  cin>>a;
    
  primo.ConfFinal_and_previous(config_f0,config_f1); 

  //For last time, operations on strings
  string _ave_ekin=locus+"ave_ekin.dat";
  char ave_ekin[_ave_ekin.size()+1];
  strcpy(ave_ekin,&_ave_ekin[0]);
  
  string _ave_epot=locus+"ave_epot.dat";
  char ave_epot[_ave_epot.size()+1];
  strcpy(ave_epot,&_ave_epot[0]);

  string _ave_etot=locus+"ave_etot.dat";
  char ave_etot[_ave_etot.size()+1];
  strcpy(ave_etot,&_ave_etot[0]);

  string _ave_temp=locus+"ave_temp.dat";
  char ave_temp[_ave_temp.size()+1];
  strcpy(ave_temp,&_ave_temp[0]);

  string _ave_press=locus+"ave_press.dat";
  char ave_press[_ave_press.size()+1];
  strcpy(ave_press,&_ave_press[0]);

  ekin.stampa(ave_ekin);
  epot.stampa(ave_epot);
  etot.stampa(ave_etot);
  temp.stampa(ave_temp);
  press.stampa(ave_press);

  /*
  string _ave_gofr=locus+"ave_gofr.dat";
  char ave_gofr[_ave_gofr.size()+1];
  strcpy(ave_gofr,&_ave_gofr[0]);

  ofstream gof_out;
  gof_out.open(ave_gofr);
  
  for(int j=0;j<nbins;j++){
    string _pesca_da="./Gofr/output_gofr_"+to_string(j)+".dat";
    char pesca_da[_pesca_da.size()+1];
    strcpy(pesca_da,&_pesca_da[0]);
    
    monteC gofr(pesca_da,primo.get_nstep()/10,blocks);//Input, elementi, numero di blocchi
    gofr.statistica();
    gof_out<<j*(primo.get_box()/2.0)/(double)nbins<<" "<<gofr.media()<<" "<<gofr.errore()<<endl;
  }
  gof_out.close();
  */
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
