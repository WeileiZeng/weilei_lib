#ifndef MM_WRITE_H
#define MM_WRITE_H

#include <itpp/itbase.h>


//int GF2mat_to_MM(GF2mat G, char* file_name="mm_temp.dat");

/** save GF2mat into MatrixMarket file*/
int GF2mat_to_MM(itpp::GF2mat G, char* file_name, int debug=0);

//int mat_to_MM(mat G, char* file_name="mm_temp.dat");

/** save (double) mat into MatrixMarket file*/
int mat_to_MM(itpp::mat G, char* file_name);



//add string compatibility

/** save GF2mat into MatrixMarket file*/
int GF2mat_to_MM(itpp::mat G, std::string file_name);

/** save (double) mat into MatrixMarket file*/
int mat_to_MM(itpp::mat G, std::string file_name);

#endif

