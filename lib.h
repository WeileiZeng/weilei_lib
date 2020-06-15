/**\file lib.h
 *\brief lib for general functions, defined within namespace common
 */

#ifndef LIB_H
#define LIB_H
#include <string> //string is already included in itpp
//#include <iostream>
//#include<fstream>
//#include <stdio.h>
#include <itpp/itbase.h>

//#include "weilei_lib.h"

/**\namespace common
 *\brief common function shared by many program
 */
namespace common{
  //general cpp util functions
  
  /** return text in red color for printing */
  std::string color_text(std::string str); 

  std::string red_text(std::string str); /**< return red std::string */

  std::string blue_text(std::string str);/**< return blue std::string */

  int get_time(int mode = 1); /**< return time in secs, milliseconds, or ... */

  //itpp based functions
  
  /** check if it is a quantum CSS code
   *@param G_x X type parity check matrix
   *@param G_z Z type parity check matrix
   */
  bool is_quantum_code(itpp::GF2mat & G_x, itpp::GF2mat & G_z);
  //check if the code is quantum or not   
  //this function is also defined elsewhere

  //moved from product.h
  /** check if it is a quantum CSS code
   *@param G_x X type parity check matrix
   *@param G_z Z type parity check matrix
   *@param C_x X type codeword generating matrix
   *@param C_z Z type codeword generating matrix
   */
  bool is_quantum_code(itpp::GF2mat &Gx,itpp::GF2mat &Gz, itpp::GF2mat &Cx,itpp::GF2mat &Cz);

  /** reduce a fat matrix with degenerate rows to a thin matrix with full rank; remove the dependent rows 
   *@param fat a GF2mat not necessary full rank
   *@return thin: a full rank GF2mat matrix
   */
  itpp::GF2mat make_it_full_rank(itpp::GF2mat fat);
  

  //int itpp::GF2matPrint(itpp::GF2mat &G, char * name = (char *) " ");
  /**print brief information of G: name, size, density
   *@param G
   *@param name name of the matrix
   */
  int GF2matPrint(itpp::GF2mat &G, std::string name);

  /**print brief information of G: name, size, density
   *@param G
   *@param name name of the matrix
   */
  int matPrint(itpp::mat G, char * name = (char *) " ");//print brief infomation of G

  /** 
   *Kronecker product of A and B 
   *@param A the outer matrix
   *@param B the inner matrix
   */
  itpp::GF2mat kron(itpp::GF2mat A, itpp::GF2mat B);//not sure how this kron is defined, maybe inversed

  /** convert int to string
   *\warning not used anywhere, should be replaced by build-in function to_string()
   */
  std::string NumberToString(int pNumber);//Convert number to std::string. //because the build in function is not supported in C++ 11

  /**< append row vectors to an itpp::GF2mat 
   *\warning This is super ineffcient, and should not be used this way
   */
  itpp::GF2mat append_vector(itpp::GF2mat G, itpp::bvec b);

  /** getitpp::GF2mat from filename with given pattern
   *full name = filename_prefix+filename_suffix
   */
  itpp::GF2mat get_GF2mat(char * filename_prefix, char * filename_suffix);

  /** getitpp::GF2mat from filename with given pattern
   *full name = parent_folder/folder/filename
   */
  itpp::GF2mat get_GF2mat(char * parent_folder, char * folder, char * filename);

  /** return density of a matrix whose first row is zero
   *remove first row and return density of the remained submatrix    
   */
  double get_error_density(itpp::GF2mat E);


  /** save mat into gnuplot data file, with a custom comment as header*/
  int mat2gnudata(itpp::mat data, std::string filename, std::string header);


}

#endif
