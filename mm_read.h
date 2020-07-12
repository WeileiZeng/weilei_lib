#ifndef MM_READ_H
#define MM_READ_H

#include <itpp/itbase.h>
//#include <string> // string is included in itbase.h


/** read MatrixMarket file and return GF2mat */
itpp::GF2mat MM_to_GF2mat(char *  file_name);
/** read MatrixMarket file and return GF2mat */
itpp::GF2mat MM_to_GF2mat(std::string  file_name);

/** read MatrixMarket file and return (double) mat */
itpp::mat MM_to_mat(char *  file_name);
/** read MatrixMarket file and return (double) mat */
itpp::mat MM_to_mat(std::string  file_name);


/** read MatrixMarket file for dense matrix and return (double) mat
 * the coordinate format is efficient for sparse matrix. for dense
 * matrix, it is saved in collumn fashion, which can be read by the
 * following file. It was included and checked in the MM_to_mat()
 * function.
 */
itpp::mat dense_MM_to_mat(char * file_name);

#endif
