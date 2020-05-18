#ifndef DIST_H
#define DIST_H
//Weilei Zeng, April 28

//#include "weilei_lib/my_lib.h"
#include <itpp/itbase.h>
#include <itpp/itcomm.h>
#include <stdio.h>
#include "my_lib.h"
#include <string>
using namespace itpp;
using namespace std;

int save_dist(int d, char * filename);

int rand_dist(GF2mat C, int perm_try=10);
//use random window decoder to find min wt of rows in C

int classical_dist(GF2mat G);
//return distance of classical binary code GH^T=0;
//G is the parity check matrix

//return H such that GH^T = 0, and rank G + rank H = n = full rank
GF2mat nullSpace(GF2mat G);

GF2mat getC(GF2mat G_x,GF2mat G_z,int flip=0);//return C_x or C_z if flip=1

int quantum_dist_v2(GF2mat G_x, GF2mat G_z, int flip=0);//without expectation value

int quantum_dist(GF2mat G_x, GF2mat G_z, int dist_expected, int flip=0);
  //right or x  distance of (G_x,G_z)
  //flip left and right if flip = 1  

int hypergraph_dist(GF2mat Aj, GF2mat Ajplus,int dist_expected,int flip=0);
//return left distance of CSS code (Aj,Ajplus^T)
//flip left and right if flip=1


//following functions for bp_new, Belief Propagation decoder.



int draw_toric_x_error(bvec error_bits);
//draw the lattice, with error bond in red

int draw_toric_x_error(bvec error_bits, string header);
//draw the lattice, with error bond in red

bvec find_error(bvec e_in, GF2mat H);
  //input: original error and parity check matrix
  //output: an error with same syndrome
  //for principle, see random window decoder

GF2mat get_check_code734(int L);//L=7n    
//return check matrix code code [7,3,4], find definition in research note.pdf    
GF2mat get_check_code743(int L);//L=7n
//return check matrix code code [7,4,3], find definition in research note.pdf
GF2mat get_check_rept(int L);//return circulant check matrix for repetition code of length L
GF2mat get_check(int generator_flag, int L);
    //return check matrix
/*switch(generator_flag){
  case 1: return get_check_rept(L);break;
  case 2: return get_check_code734(L);break;
  case 3: return get_check_code743(L);break;*/

LDPC_Code GF2mat_to_LDPC_Code(GF2mat G);

LDPC_Code MM_to_LDPC_Code(char * filename);
  //convert GF2mat saved in .mm file to LDPC_Code  
    

#endif

