/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <iostream>     // cin, cout: Standard Input/Output Streams Library
#include <fstream>      // Stream class to both read and write from/to files.
#include <cmath>        // rint, pow
#include <string>      
#include "random.h"     // generating random numbers
#include "MolDyn_NVE.h"

using namespace std;

MolDyn::MolDyn(const char* input, const char* config){ //Prepare all stuff for the simulation

/////////////////////////////////////////////////////////////////////////
//Initialize random number generators
   Random rnd;
   int seed[4];
   int p1, p2;
   ifstream Primes("Primes");
   if (Primes.is_open()){
      Primes >> p1 >> p2 ;
   } else cerr << "PROBLEM: Unable to open Primes" << endl;
   Primes.close();

   ifstream read_input("seed.in");
   string property;
   if (read_input.is_open()){
      while ( !read_input.eof() ){
         read_input >> property;
         if( property == "RANDOMSEED" ){
            read_input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
            rnd.SetRandom(seed,p1,p2);
         }
      }
      read_input.close();
   } else cerr << "PROBLEM: Unable to open seed.in" << endl;
/////////////////////////////////////////////////////////////////////////

  ifstream ReadInput,ReadConf;
  //double ep, ek, pr, et, vir;

  cout << "Classic Lennard-Jones fluid        " << endl;
  cout << "Molecular dynamics simulation in NVE ensemble  " << endl << endl;
  cout << "Interatomic potential v(r) = 4 * [(1/r)^12 - (1/r)^6]" << endl << endl;
  cout << "The program uses Lennard-Jones units " << endl;

  ReadInput.open(input); //Read input

  ReadInput >> temp;

  ReadInput >> npart;
  cout << "Number of particles = " << npart << endl;

  ReadInput >> rho;
  cout << "Density of particles = " << rho << endl;
  vol = (double)npart/rho;
  cout << "Volume of the simulation box = " << vol << endl;
  box = pow(vol,1.0/3.0);
  cout << "Edge of the simulation box = " << box << endl;

  ReadInput >> rcut;
  ReadInput >> delta;
  ReadInput >> nstep;
  ReadInput >> iprint;

  cout << "The program integrates Newton equations with the Verlet method " << endl;
  cout << "Time step = " << delta << endl;
  cout << "Number of steps = " << nstep << endl << endl;
  ReadInput.close();
  
  //Inizialize arrays
  x=new double[npart];
  y=new double[npart];
  z=new double[npart];
  xold=new double[npart];
  yold=new double[npart];
  zold=new double[npart];  
  vx=new double[npart];
  vy=new double[npart];
  vz=new double[npart];

//Read initial configuration
  cout << "Read initial configuration from file "<<config << endl << endl;
  ReadConf.open(config);
  for (int i=0; i<npart; ++i){
    ReadConf >> x[i] >> y[i] >> z[i];
    x[i] = x[i] * box;
    y[i] = y[i] * box;
    z[i] = z[i] * box;
  }
  ReadConf.close();

//Prepare initial velocities
   cout << "Prepare random velocities with center of mass velocity equal to zero " << endl << endl;
   double sumv[3] = {0.0, 0.0, 0.0};
   for (int i=0; i<npart; ++i){
     vx[i] = rnd.Gauss(0,sqrt(temp));//Da Maxwell-Boltzmann in unita' di Lennard-Jones
     vy[i] = rnd.Gauss(0,sqrt(temp));
     vz[i] = rnd.Gauss(0,sqrt(temp));

     sumv[0] += vx[i];
     sumv[1] += vy[i];
     sumv[2] += vz[i];
   }
   for (int idim=0; idim<3; ++idim) sumv[idim] /= (double)npart;
   double sumv2 = 0.0, fs;
   for (int i=0; i<npart; ++i){
     vx[i] = vx[i] - sumv[0];
     vy[i] = vy[i] - sumv[1];
     vz[i] = vz[i] - sumv[2];

     sumv2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
   }
   sumv2 /= (double)npart;

   fs = sqrt(3 * temp / sumv2);   // fs = velocity scale factor

   cout<<"Scale factor = " <<fs<<endl<<endl;
   
   for (int i=0; i<npart; ++i){
     vx[i] *= fs;
     vy[i] *= fs;
     vz[i] *= fs;

     xold[i] = Pbc(x[i] - vx[i] * delta);
     yold[i] = Pbc(y[i] - vy[i] * delta);
     zold[i] = Pbc(z[i] - vz[i] * delta);
     /////////////////////////////////////////////////
     //Ricontrolla!!!
   }

   //bin_size=(box/2.0)/(double)nbins;bin_size=(box/2.0)/(double)nbins;
}

MolDyn::MolDyn(const char* input, const char* config1, const char* config2){ //Prepare all stuff for the simulation

  ifstream ReadInput,ReadConf;
  //double ep, ek, pr, et, vir;

  cout << "Classic Lennard-Jones fluid        " << endl;
  cout << "Molecular dynamics simulation in NVE ensemble  " << endl << endl;
  cout << "Interatomic potential v(r) = 4 * [(1/r)^12 - (1/r)^6]" << endl << endl;
  cout << "The program uses Lennard-Jones units " << endl;

  ReadInput.open(input); //Read input

  ReadInput >> temp;

  ReadInput >> npart;
  cout << "Number of particles = " << npart << endl;

  ReadInput >> rho;
  cout << "Density of particles = " << rho << endl;
  vol = (double)npart/rho;
  cout << "Volume of the simulation box = " << vol << endl;
  box = pow(vol,1.0/3.0);
  cout << "Edge of the simulation box = " << box << endl;

  ReadInput >> rcut;
  ReadInput >> delta;
  ReadInput >> nstep;
  ReadInput >> iprint;

  cout << "The program integrates Newton equations with the Verlet method " << endl;
  cout << "Time step = " << delta << endl;
  cout << "Number of steps = " << nstep << endl << endl;
  ReadInput.close();
  
  //Inizialize arrays
  x=new double[npart];
  y=new double[npart];
  z=new double[npart];
  xold=new double[npart];
  yold=new double[npart];
  zold=new double[npart];  
  vx=new double[npart];
  vy=new double[npart];
  vz=new double[npart];

//Read initial configuration
  cout << "Read initial configuration from file "<<config1 << endl << endl;
  ReadConf.open(config1);
  for (int i=0; i<npart; ++i){
    ReadConf >> xold[i] >> yold[i] >> zold[i];
    xold[i] = xold[i] * box;
    yold[i] = yold[i] * box;
    zold[i] = zold[i] * box;
  }
  ReadConf.close();

//Read next-step configuration
  cout << "Read successive configuration from file "<<config2 << endl << endl;
  ReadConf.open(config2);
  for (int i=0; i<npart; ++i){
    ReadConf >> x[i] >> y[i] >> z[i];
    x[i] = x[i] * box;
    y[i] = y[i] * box;
    z[i] = z[i] * box;
  }
  ReadConf.close();
  ////////////////////////////////////////////////////////////
//Calculate velocities
  cout << "Calculate velocities for the inserted configuration of particles" << endl << endl;
  
  //Per prima cosa faccio fare provvisoriamente un passo al sistema
  double xnew, ynew, znew;
  double* fx=new double[npart];
  double* fy=new double[npart];
  double* fz=new double[npart];

  double sumv[3] = {0.0, 0.0, 0.0};//Uso questi vettori per calcolare la temperatura partendo dalle velocità
    
  for(int i=0; i<npart; ++i){ //Force acting on particle i
    fx[i] = Force(i,0);
    fy[i] = Force(i,1);
    fz[i] = Force(i,2);
  }

  for(int i=0; i<npart; ++i){ //Verlet integration scheme

    xnew = Pbc( 2.0 * x[i] - xold[i] + fx[i] * pow(delta,2) );
    ynew = Pbc( 2.0 * y[i] - yold[i] + fy[i] * pow(delta,2) );
    znew = Pbc( 2.0 * z[i] - zold[i] + fz[i] * pow(delta,2) );

    vx[i] = Pbc(xnew - xold[i])/(2.0 * delta);
    vy[i] = Pbc(ynew - yold[i])/(2.0 * delta);
    vz[i] = Pbc(znew - zold[i])/(2.0 * delta);

    sumv[0] += vx[i];
    sumv[1] += vy[i];
    sumv[2] += vz[i];
  }
   
   for (int idim=0; idim<3; ++idim) sumv[idim] /= (double)npart;
   double sumv2 = 0.0, fs;
   for (int i=0; i<npart; ++i){
     vx[i] = vx[i] - sumv[0];
     vy[i] = vy[i] - sumv[1];
     vz[i] = vz[i] - sumv[2];

     sumv2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
   }
   sumv2 /= (double)npart;

   fs = sqrt(3 * temp / sumv2);   // fs = velocity scale factor

   if(abs(fs-1)>0.04){
   cout<<"Scale factor = " <<fs<<endl<<endl;
   
   for (int i=0; i<npart; ++i){
     vx[i] *= fs;
     vy[i] *= fs;
     vz[i] *= fs;

     xold[i] = Pbc(x[i] - vx[i] * delta);
     yold[i] = Pbc(y[i] - vy[i] * delta);
     zold[i] = Pbc(z[i] - vz[i] * delta);
   }
   }
   else{
     cout<<"No need for scale factor"<<endl<<endl;
   }
   //bin_size=(box/2.0)/(double)nbins;
}

void MolDyn::Move(void){ //Move particles with Verlet algorithm
  double xnew, ynew, znew;
  double* fx=new double[npart];
  double* fy=new double[npart];
  double* fz=new double[npart];

  for(int i=0; i<npart; ++i){ //Force acting on particle i
    fx[i] = Force(i,0);
    fy[i] = Force(i,1);
    fz[i] = Force(i,2);
  }

  for(int i=0; i<npart; ++i){ //Verlet integration scheme

    xnew = Pbc( 2.0 * x[i] - xold[i] + fx[i] * pow(delta,2) );
    ynew = Pbc( 2.0 * y[i] - yold[i] + fy[i] * pow(delta,2) );
    znew = Pbc( 2.0 * z[i] - zold[i] + fz[i] * pow(delta,2) );

    vx[i] = Pbc(xnew - xold[i])/(2.0 * delta);
    vy[i] = Pbc(ynew - yold[i])/(2.0 * delta);
    vz[i] = Pbc(znew - zold[i])/(2.0 * delta);

    xold[i] = x[i];
    yold[i] = y[i];
    zold[i] = z[i];

    x[i] = xnew;
    y[i] = ynew;
    z[i] = znew;
  }
  return;
}

double MolDyn::Force(int ip, int idir){ //Compute forces as -Grad_ip V(r)
  double f=0.0;
  double dvec[3], dr;

  for (int i=0; i<npart; ++i){
    if(i != ip){
      dvec[0] = Pbc( x[ip] - x[i] );  // distance ip-i in pbc
      dvec[1] = Pbc( y[ip] - y[i] );
      dvec[2] = Pbc( z[ip] - z[i] );

      dr = dvec[0]*dvec[0] + dvec[1]*dvec[1] + dvec[2]*dvec[2];
      dr = sqrt(dr);

      if(dr < rcut){
        f += dvec[idir] * (48.0/pow(dr,14) - 24.0/pow(dr,8)); // -Grad_ip V(r)
      }
    }
  }
  
  return f;
}

void MolDyn::Measure(){ //Properties measurement
  //int bin;
  double v, K, vij, virial_ij, W;
  double dx, dy, dz, dr;
  ofstream Epot, Ekin, Etot, Temp, Press;//, Gofr;

  Epot.open("output_epot.dat",ios::app);
  Ekin.open("output_ekin.dat",ios::app);
  Temp.open("output_temp.dat",ios::app);
  Etot.open("output_etot.dat",ios::app);
  Press.open("output_press.dat",ios::app);

  v = 0.0; //reset observables
  K = 0.0;
  W=0.0;

  /*
  for(int j=0;j<nbins;j++){
    bins[j]=0;
  }
  */
  
  //cycle over pairs of particles
  for (int i=0; i<npart-1; ++i){
    for (int j=i+1; j<npart; ++j){
      
      dx = Pbc( x[i] - x[j] );
      dy = Pbc( y[i] - y[j] );
      dz = Pbc( z[i] - z[j] );
      
      dr = dx*dx + dy*dy + dz*dz;
      dr = sqrt(dr);

      /*
      //update of the histogram of g(r)
      if(dr<box/2){
	//cout<<int((dr/(box/2))*nbins)+igofr<<endl;
	bins[int(dr/bin_size)]=bins[int(dr/bin_size)]+2;
      }
      */
     
      if(dr < rcut){
	vij = 4.0/pow(dr,12) - 4.0/pow(dr,6);
	virial_ij= 48.0/pow(dr,12) - 24.0/pow(dr,6);

	//Potential energy
	v += vij;
	W+=virial_ij;
      }
    }          
  }

  //Term for pressure
  W/=npart;

//Kinetic energy
  for (int i=0; i<npart; ++i) K += 0.5 * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);
   
    stima_pot = v/(double)npart; //Potential energy
    stima_kin = K/(double)npart; //Kinetic energy
    stima_temp = (2.0 / 3.0) * K/(double)npart; //Temperature
    stima_etot = (K+v)/(double)npart; //Total enery
    stima_press = rho*temp+W/(3*vol); //Pressure

    Epot << stima_pot  << endl;
    Ekin << stima_kin  << endl;
    Temp << stima_temp << endl;
    Etot << stima_etot << endl;
    Press << stima_press << endl;

    /*
    for(int j=0;j<nbins;j++){
      Gofr.open("./Gofr/output_gofr_"+to_string(j)+".dat",ios::app);
      double V_r_dr=4*M_PI/3*( pow((j+1)*bin_size,3)-pow((j)*bin_size,3) );
      Gofr << double(bins[j])/(npart*rho*V_r_dr)<<endl;
      Gofr.close();
    }
    */

    Epot.close();
    Ekin.close();
    Temp.close();
    Etot.close();
    Press.close();
    
    return;
}


void MolDyn::ConfFinal(const char* output){ //Write final configuration
  ofstream WriteConf;

  cout <<endl<< "Print final configuration to file "<< output << endl << endl;
  WriteConf.open(output);

  for (int i=0; i<npart; ++i){
    WriteConf << x[i]/box << "   " <<  y[i]/box << "   " << z[i]/box << endl;
  }
  WriteConf.close();
  return;
}

void MolDyn::ConfFinal_and_previous(const char* output1, const char* output2){ //Write final configuration, make the system evolve once more and write it a second time

  ConfFinal(output1);//Print output
  Move();//Make the system evolve
  cout << "Printing also the next step." << endl;
  ConfFinal(output2);//Print output
  
  return;
}

void MolDyn::ConfXYZ(int nconf){ //Write configuration in .xyz format
  ofstream WriteXYZ;

  WriteXYZ.open("frames/config_" + to_string(nconf) + ".xyz");
  WriteXYZ << npart << endl;
  WriteXYZ << "This is only a comment!" << endl;
  for (int i=0; i<npart; ++i){
    WriteXYZ << "LJ  " << Pbc(x[i]) << "   " <<  Pbc(y[i]) << "   " << Pbc(z[i]) << endl;
  }
  WriteXYZ.close();
}

double MolDyn::Pbc(double r){  //Algorithm for periodic boundary conditions with side L=box
    return r - box * rint(r/box);
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
