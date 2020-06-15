#ifndef DIST_H
#define DIST_H
//Weilei Zeng, April 28

//#include "weilei_lib/my_lib.h"
//#include "weilei_lib.h"


//#include <string> //already inlcuded in itbase.h
//#include <itpp/itbase.h> //already included in itcomm.h
#include <itpp/itcomm.h> //LDPC_Code
//#include <stdio.h>

namespace common {

//const int MAX_M=6;//maximum of the length of the complex chain
const int INF=999;//infinity distance


int min_wt_decoding(itpp::GF2mat C);

//G for gauge operators, and C for bare logical operators
//code word c = alpha_C*C+alpha_G*G, where alpha_C \neq 0
int min_wt_decoding(itpp::GF2mat C, itpp::GF2mat G);

int save_dist(int d, char * filename);

int rand_dist(itpp::GF2mat C, int perm_try=10);
//use random window decoder to find min wt of rows in C

int classical_dist(itpp::GF2mat G);
//return distance of classical binary code GH^T=0;
//G is the parity check matrix

//return H such that GH^T = 0, and rank G + rank H = n = full rank
itpp::GF2mat nullSpace(itpp::GF2mat G);

itpp::GF2mat getC(itpp::GF2mat G_x,itpp::GF2mat G_z,int flip=0);//return C_x or C_z if flip=1

int quantum_dist_v2(itpp::GF2mat G_x, itpp::GF2mat G_z, int flip=0);//without expectation value

int quantum_dist(itpp::GF2mat G_x, itpp::GF2mat G_z, int dist_expected, int debug, int flip=0);
  //right or x  distance of (G_x,G_z)
  //flip left and right if flip = 1  

int hypergraph_dist(itpp::GF2mat Aj, itpp::GF2mat Ajplus,int dist_expected,int flip=0);
//return left distance of CSS code (Aj,Ajplus^T)
//flip left and right if flip=1


//following functions for bp_new, Belief Propagation decoder.



int draw_toric_x_error(itpp::bvec error_bits);
//draw the lattice, with error bond in red

int draw_toric_x_error(itpp::bvec error_bits, std::string header);
//draw the lattice, with error bond in red

itpp::bvec find_error(itpp::bvec e_in, itpp::GF2mat H);
  //input: original error and parity check matrix
  //output: an error with same syndrome
  //for principle, see random window decoder

itpp::GF2mat get_check_code734(int L);//L=7n    
//return check matrix code code [7,3,4], find definition in research note.pdf    
itpp::GF2mat get_check_code743(int L);//L=7n
//return check matrix code code [7,4,3], find definition in research note.pdf
itpp::GF2mat get_check_rept(int L);//return circulant check matrix for repetition code of length L
itpp::GF2mat get_check(int generator_flag, int L);
    //return check matrix
/*switch(generator_flag){
  case 1: return get_check_rept(L);break;
  case 2: return get_check_code734(L);break;
  case 3: return get_check_code743(L);break;*/

itpp::LDPC_Code GF2mat_to_LDPC_Code(itpp::GF2mat G);

itpp::LDPC_Code MM_to_LDPC_Code(char * filename);
  //convert GF2mat saved in .mm file to LDPC_Code  
    

} //namespace common
#endif

