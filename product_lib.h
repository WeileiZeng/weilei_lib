//#ifndef CONCATENATION_LIB_H
//#define CONCATENATION_LIB_H
#ifndef PRODUCT_LIB_H
#define PRODUCT_LIB_H

#include "dist.h"
#include <itpp/itbase.h>
//#include <itpp/itcomm.h>
//#include <stdio.h>
#include "weilei_lib.h"


const int MAX_M=6;//maximum of the length of the complex chain
//const int INF=999;//infinity distance




int getRandomQuantumCode(int n,int Gx_row,int Gz_row, itpp::GF2mat &Gx,itpp::GF2mat &Gz, itpp::GF2mat &Cx,itpp::GF2mat &Cz);

int getGoodQuantumCode(int n,int Gx_row,int Gz_row, itpp::GF2mat &Gx,itpp::GF2mat &Gz, itpp::GF2mat &Cx,itpp::GF2mat &Cz, int debug);


void set_submatrix(itpp::GF2mat & G, itpp::GF2mat sub, int row, int col);

int generate_code(itpp::GF2mat & Gax, itpp::GF2mat & Gaz, int na, int Gax_row, int id_Gax, int Gaz_row, int id_Gaz, int debug);

//files used in concatenated codes and product codes.

//int reduce(itpp::GF2mat Gax, itpp::GF2mat Gaz, itpp::GF2mat Gbx, itpp::GF2mat Gbz,int ddax,int ddaz,int ddbx,int ddbz);
//int concatenate(itpp::GF2mat Gax, itpp::GF2mat Gaz, itpp::GF2mat Gbx, itpp::GF2mat Gbz,int ddax,int ddaz,int ddbx,int ddbz);



// a version include both reduce and concatenation
// mode=1 for reduce/subsystem product
// mode=2 for concatenation
//only dz is checked cause dx is known to be tight
int product(itpp::GF2mat Gax, itpp::GF2mat Gaz, itpp::GF2mat Gbx, itpp::GF2mat Gbz,int ddax,int ddaz,int ddbx,int ddbz, int debug, int mode);




#endif //PRODUCT_LIB_H
