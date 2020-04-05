#ifndef DIST_H
#define DIST_H
//Weilei Zeng, April 28

//#include "weilei_lib/my_lib.h"
#include <itpp/itbase.h>
#include <itpp/itcomm.h>
#include <stdio.h>
#include "my_lib.h"
using namespace itpp;
using namespace std;

int save_dist(int d, char * filename);

int rand_dist(GF2mat C, int perm_try=10);
//use random window decoder to find min wt of rows in C

int classical_dist(GF2mat G);
//return distance of classical binary code GH^T=0;
//G is the parity check matrix

GF2mat getC(GF2mat G_x,GF2mat G_z,int flip=0);//return C_x or C_z if flip=1

int quantum_dist_v2(GF2mat G_x, GF2mat G_z, int flip=0);//without expectation value

int quantum_dist(GF2mat G_x, GF2mat G_z, int dist_expected, int flip=0);
  //right or x  distance of (G_x,G_z)
  //flip left and right if flip = 1  

int hypergraph_dist(GF2mat Aj, GF2mat Ajplus,int dist_expected,int flip=0);
//return left distance of CSS code (Aj,Ajplus^T)
//flip left and right if flip=1

#endif

