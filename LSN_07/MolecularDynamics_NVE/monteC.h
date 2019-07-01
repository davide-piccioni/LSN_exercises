/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#ifndef __monteC__
#define __monteC__

#include "random.h"

class monteC {
  
 private:
  unsigned int N; //Numero di lanci totali
  unsigned int k;//Numero di blocchi
  unsigned int L; //Elemeni per blocco, sarà L=N/k
  double * raw;//Vettore che contiene i valori simulati
  double * prog_1;//Vettore delle medie
  double * prog_2;//Vettore delle medie dei quadrati
  double * err_prog;//Vettore delle incertezze
  Random rnd;//Generatore
  
protected:

public:
  // costruttore di default: N=100000,k=100
  monteC();
  // costruttore: N,k
  monteC(unsigned int N_tuo, unsigned int k_tuo);
  // costruttore: prende gia' un vettore con raw riempito e ne copia gli elementi
  monteC(double * my_raw, unsigned int N_tuo, unsigned int k_tuo);
  // costruttore da file (devi gia' sapere il numero di elementi da inserire)
  monteC(char * input_file, unsigned int N_tuo, unsigned int k_tuo);

  // destructor
  ~monteC();
  // methods
  unsigned int get_N()const{return N;}
  unsigned int get_k()const{return k;}

  void set_raw(double * my_raw);//Inserisci tu il vettore raw

  double**  statistica();//Calcola k valori medi e i rispettivi errori (uno per ciascun blocco), restituendoli in una matrice 2 x N
  void stampa_schermo()const; //Stampa a video: prima colonna valori, seconda errori
  void stampa(const char* output)const; //Stampa su file: prima colonna valori, seconda errori
  double media(){return prog_1[k-1];}//Rende la media piu' precisa
  double errore(){return err_prog[k-1];}//Rende l'errore minimo
 
};

#endif // __monteC__

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
