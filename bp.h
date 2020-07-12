#ifndef BP_H
#define BP_H

#include <itpp/itcomm.h> //qllr_to_bvec


itpp::bvec reduce_weight(itpp::bvec e, itpp::GF2mat G);
//reduce the weight of en error vector e, G include all cycles (equivalent error). G is assument to be the generator matrix of   an Quantum LDPC CSS code. such that rows of G correspond to the smallest cycles. It is assument that e is an effectively small/zero weight error, which contains many trivial cycles. This function will reduce most of those small individual cycles so that we can count the effective weight of e.
//  complexity of this program: ~ G.rows() * weight(e)   


itpp::bvec qllr_to_bvec(itpp::QLLRvec llr, int bound);
  //find the appropriate bound value and return binary vector from this qllr.
  //Its sign is the same as the sign of "bound"    

int check_matrices(itpp::GF2mat * G, itpp::GF2mat * H, itpp::GF2mat *U, itpp::GF2mat * W, itpp::mat * K);
  //check if those matrices are valid 

//for partial sum decoder
//read those matrices in the folder
int read_matrices_for_partial_sum(char * filename_prefix, char * filename_suffix, itpp::GF2mat * G, itpp::GF2mat * H, itpp::GF2mat *U, itpp::GF2mat * W, itpp::mat * K);

int remove_rows(itpp::GF2mat *G, itpp::bvec rows_to_remove);

#endif

