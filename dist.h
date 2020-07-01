/** \file dist.h
 *\brief distance related functions, defined within namespace common
 */

#ifndef DIST_H
#define DIST_H
//Weilei Zeng, April 28


//#include <string> //already inlcuded in itbase.h
//#include <itpp/itbase.h> //already included in itcomm.h
#include <itpp/itcomm.h> //LDPC_Code
//#include <stdio.h>

namespace common {

  //const int MAX_M=6;//maximum of the length of the complex chain

  const int INF=999;/**< Infinite distance in random window decoder*/

  /** min weight decoder for classical code
   * @param C codeword generating matrix
   * @return min weight of all codewords (any combination of rows in the generating matrix)
   */
  int min_wt_decoding(itpp::GF2mat C);

  /** min weight decoder for CSS code
   * distance The distance between \f$ y=ax^2+b \f$ 
   * @param C codeword generating matrix
   * @param G gauge matrix
   * @return \f$ d_{\text{min}} = \text{min}_{\alpha \neq 0, \beta} \text{wgt} (\alpha C+ \beta G), \alpha \neq 0 \f$ 
   */
  int min_wt_decoding(itpp::GF2mat C, itpp::GF2mat G);

  /** save distance to a file, as a single-element itpp::mat */
  int save_dist(int d, char * filename);

  /** use random window decoder to find min wt of rows in C.
   * min_wt_decoding is used when C.rows()<7
   *@param C codeword generating matrix
   *@param perm_try=10 number of random trials
   */
  int rand_dist(itpp::GF2mat C, int perm_try=10);

  /** use random window decoder to find distance of a classical code with parity check patrix G
   *@param G parity check matrix \f$GC^T=0\f$
   * min_wt_decoding is used when C.rows()<7
   * perm_try=10 always. No interface is given to change it yet.
   */
  int classical_dist(itpp::GF2mat G);


  /** @return H such that \f$GH^T = 0\f$, and rank G + rank H = n = full rank */
  itpp::GF2mat nullSpace(itpp::GF2mat G);

  /** get codeword generating matrix for CSS codes
   *@param G_x X type party check matrix
   *@param G_z Z type party check matrix
   *@param flip=0 whether flipping or nor
   *@return C_x when flip=0
   *@return C_z when flip=1
   */
  itpp::GF2mat getC(itpp::GF2mat G_x,itpp::GF2mat G_z,int flip=0);//return C_x or C_z if flip=1

  /** get estimated distance for CSS codes
   *@param G_x X type party check matrix
   *@param G_z Z type party check matrix
   *@param flip=0 whether flipping or nor
   *@return d_x when flip=0
   *@return d_z when flip=1
   */
  int quantum_dist_v2(itpp::GF2mat G_x, itpp::GF2mat G_z, int flip=0);//without expectation value

    /** get estimated distance for CSS codes
   *@param G_x X type party check matrix
   *@param G_z Z type party check matrix
   *@param flip=0 whether flipping or nor
   *@param dist_expected use expected distance to control number of trials
   *@return d_x when flip=0
   *@return d_z when flip=1
   */
  int quantum_dist(itpp::GF2mat G_x, itpp::GF2mat G_z, int dist_expected, int debug, int flip=0);
  //right or x  distance of (G_x,G_z)
  //flip left and right if flip = 1  

  /** get estimated distance of a CSS code generated from a chain complex \f$ A_jA_{j+1}=0 \f$
   *@param Aj
   *@param Ajplus
   *@param dist_expected use expected distance to control number of trials
   *@param flip=0 (default) or 1, flip left and right if flip=1
   *@return left distance of CSS code (Aj,Ajplus^T) 
   */
  int hypergraph_dist(itpp::GF2mat Aj, itpp::GF2mat Ajplus,int dist_expected,int flip=0);


  //following functions for bp_new, Belief Propagation decoder.


  /** draw the lattice, with error bond in red */
  int draw_toric_x_error(itpp::bvec error_bits);

  /** a wrapper of draw_toric_x_error, with an extra header*/
  int draw_toric_x_error(itpp::bvec error_bits, std::string header);

  /** find error with same syndrome for a classical code, can be used for CSS codes as well
   *@param e_in original error
   *@param H  parity check matrix
   *@return an error with same syndrome (It is not necessarily an equivalent error)
  //for principle, see random window decoder
  */
  itpp::bvec find_error(itpp::bvec e_in, itpp::GF2mat H);

  /** return check matrix code code [7,3,4], find definition in research note.pdf
   *@param L size of the code, must be a multiple of 7
   */
  itpp::GF2mat get_check_code734(int L);//L=7n    

  /** return check matrix code code [7,4,3], find definition in research note.pdf 
   *@param L size of the code, must be a multiple of 7*/
  itpp::GF2mat get_check_code743(int L);//L=7n
  /** return circulant check matrix for repetition code of length L */
  itpp::GF2mat get_check_rept(int L);

  /** return check matrix
   *@param L size of the code
   switch(generator_flag){
   case 1: return get_check_rept(L);break;
   case 2: return get_check_code734(L);break;
   case 3: return get_check_code743(L);break;
  */  
  itpp::GF2mat get_check(int generator_flag, int L);


  itpp::LDPC_Code GF2mat_to_LDPC_Code(itpp::GF2mat G);

  /** convert GF2mat saved in .mm file to LDPC_Code  */
  itpp::LDPC_Code MM_to_LDPC_Code(char * filename);
  
    

} //namespace common

#endif //#ifndef DIST_H

