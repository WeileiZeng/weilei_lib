#ifndef LIB_H
#define LIB_H
#include <string>
#include <itpp/itbase.h>
using namespace itpp;
using namespace std;

bool is_quantum_code(GF2mat G_x, GF2mat G_z);
  //check if the code is quantum or not   
GF2mat make_it_full_rank(GF2mat fat);
  //reduce a fat matrix with degenerate rows to a thin matrix with full rank; remove the dependent rows 
int GF2matPrint(GF2mat G, char * name);
GF2mat kron(GF2mat A, GF2mat B);
string NumberToString(int pNumber);
GF2mat append_vector(GF2mat G,bvec b);
GF2mat get_GF2mat(char * filename_prefix, char * filename_suffix);
GF2mat get_GF2mat(char * parent_folder, char * folder, char * filename);
double get_error_density(GF2mat E);
#endif
