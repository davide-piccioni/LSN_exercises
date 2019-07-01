/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#ifndef __TSP_annealing__
#define __TSP_annealing__

#include "random.h"
#include <vector>

using namespace std;

struct pathway{
  vector<int> alleles;
  double L;
};

class anneal_tsp {
  
 private:
  int N; //Numero di citta'
  pathway our_path; //Il nostro percorso
  vector<double> cities_x;
  vector<double> cities_y; 
  Random rnd; //Generatore
  void calculate_L();

 public:
  // Costruttore da un input_file che contiene come primo elemento il numero di citta', poi le posizioni delle citta' nel piano
  anneal_tsp(const char* file_input);

  // Costruttore analogo a quello sopra, che come secondo valore prende la riga da leggere nel file "Primes" per inizializzare il generatore di numeri casuali. Controllare quante righe ha il file "Primes".
  anneal_tsp(const char* file_input,int n_riga);

  // Distruttore
  ~anneal_tsp();

  //Gives 0 if everything is ok, 1 otherwise
  bool check_path();

  // Returns L
  double get_L();

  // Print the current pathway.
  void print_path_numbers();

  //Print the coordinates of the cities in the order the pathway tells you to visit them, print in the file you give
  void print_route(const char*);
  
  
  //MUTATION METHODS:

  //They will change the pathway according to Metropolis algorithm with temperature beta and L(x) as energy.
  //They return 0 if the mutation hasn't happened, 1 otherwise
  //Gets single exchanges
  bool exchange(double beta);
  //Local shift
  bool local_shift(double beta);
  //Local permutation
  bool local_permut(double beta);


  //ACHTUNG: This will change the pathway for n_steps times according to Metropolis algorithm with temperature beta and L(x) as energy.
  //The moves are proposed randomly between exchange, local_shift and local_permutation
  //It returns a vector of 3 elements with the acceptance rate respectevly of exchange, local_shift and local_permutation
  vector<double> move(double beta, int n_steps);
  
 
};

#endif // __TSP_annealing__

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
