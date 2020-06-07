#ifndef MM_READ_H
#define MM_READ_H

#include <itpp/itbase.h>
#include <string>
using namespace std;
using namespace itpp;


GF2mat MM_to_GF2mat(char *  file_name);
GF2mat MM_to_GF2mat(string  file_name);


mat MM_to_mat(char *  file_name);
mat MM_to_mat(string  file_name);

//the coordinate format is efficient for sparse matrix. for dense matrix, it is saved in collomn fashion, which can be read by the following file. It was included and checked in the MM_to_mat() function.
mat dense_MM_to_mat(char * file_name);
#endif
