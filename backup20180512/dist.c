//Weilei Zeng, April 28
//return min weight of rows in C, which is the distance of the code
//random window method is applied with default number of permutation 10.

//#include "weilei_lib/my_lib.h"
#include "dist.h"
#include <itpp/itbase.h>
#include <itpp/itcomm.h>
#include <stdio.h>
#include "my_lib.h"
using namespace itpp;
using namespace std;
int save_dist(int d,char * filename){
  mat mat_d(1,1);
  mat_d.set(0,0,d);
  mat_to_MM(mat_d,filename);
  return 0;
}

int rand_dist(GF2mat C, int perm_try){
  //use random window decoder to find min wt of C
  RNG_randomize();
  bvec row_vec,zero=zeros_b(C.cols());
  int wt,min_wt=C.cols();
  //  int perm_try=10;//need optimization?
  ivec perm;
  GF2mat T,U;
  ivec P;
  for (int j=0;j<perm_try;j++){
    //cout<<"perm_j = "<<j<<endl;
    perm = sort_index( randu(C.cols()) );
    C.permute_cols(perm,false);
    //no need to permute back; can also permute the rows
    C.T_fact(T,U,P);
    for (int i = 0;i<C.rows();i++){
      row_vec = C.get_row(i);
      wt = BERC::count_errors(zero,row_vec);
      if (wt < min_wt){
	min_wt = wt;
	//	cout<<"got new min wt = "<<min_wt<<" when j = "<<j<<endl;
      }
    }
  }
  return min_wt;
}
