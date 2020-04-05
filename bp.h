#ifndef BP_H
#define BP_H
//Weilei Zeng, April 4, 2020


#include <itpp/itbase.h>
#include <itpp/itcomm.h>
#include <stdio.h>
#include "my_lib.h"
using namespace itpp;
using namespace std;


bvec reduce_weight(bvec e, GF2mat G);
//reduce the weight of en error vector e, G include all cycles (equivalent error). G is assument to be the generator matrix of   an Quantum LDPC CSS code. such that rows of G correspond to the smallest cycles. It is assument that e is an effectively small/zero weight error, which contains many trivial cycles. This function will reduce most of those small individual cycles so that we can count the effective weight of e.
//  complexity of this program: ~ G.rows() * weight(e)   


bvec qllr_to_bvec(QLLRvec llr, int bound);
  //find the appropriate bound value and return binary vector from this qllr.
  //Its sign is the same as the sign of "bound"    

int check_matrices(GF2mat * G, GF2mat * H, GF2mat *U, GF2mat * W, mat * K);
  //check if those matrices are valid 

//for partial sum decoder
//read those matrices in the folder
int read_matrices_for_partial_sum(char * filename_prefix, char * filename_suffix, GF2mat * G, GF2mat * H, GF2mat *U, GF2mat * W, 
				  mat * K);

#endif

