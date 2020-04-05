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



#endif

