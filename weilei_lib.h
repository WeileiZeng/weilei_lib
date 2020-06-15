/** \file weilei_lib.h
 *\brief this file links all other headfiles in this folder.
 */

#ifndef WEILEI_LIB_H
#define WEILEI_LIB_H

//no namespace
#include "mm_read.h" //file to GF2mat
#include "mm_write.h" //GF2mat to file

//namespace common
#include "lib.h" //int to string
#include "dist.h" //random window method to clculate distance

//#include "concatenation_lib.h"
#include "product_lib.h"
#include "bp.h"

//#include <math.h> //sqrt

//#include "bp_decoder.h"


//these two constant should not be here
//const int MAX_M=6;//maximum of the length of the complex chain
//const int INF=999;//infinity distance
//const double INF_BP = 1000000;
#endif
