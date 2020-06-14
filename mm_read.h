#ifndef MM_READ_H
#define MM_READ_H

#include <itpp/itbase.h>
#include <string>
//using namespace std;
//using namespace itpp;


itpp::GF2mat MM_to_GF2mat(char *  file_name);

itpp::GF2mat MM_to_GF2mat(std::string  file_name);


itpp::mat MM_to_mat(char *  file_name);
itpp::mat MM_to_mat(std::string  file_name);

//the coordinate format is efficient for sparse matrix. for dense matrix, it is saved in collomn fashion, which can be read by the following file. It was included and checked in the MM_to_mat() function.
itpp::mat dense_MM_to_mat(char * file_name);
#endif
