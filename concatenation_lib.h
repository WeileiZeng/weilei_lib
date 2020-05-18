#ifndef CONCATENATION_LIB_H
#define CONCATENATION_LIB_H
#include "dist.h"
#include <itpp/itbase.h>
#include <itpp/itcomm.h>
#include <stdio.h>
#include "my_lib.h"
using namespace itpp;
using namespace std;


// same functionis also defined in lib.h
bool is_quantum_code(GF2mat &Gx,GF2mat &Gz, GF2mat &Cx,GF2mat &Cz);


int getRandomQuantumCode(int n,int Gx_row,int Gz_row, GF2mat &Gx,GF2mat &Gz, GF2mat &Cx,GF2mat &Cz);

int getGoodQuantumCode(int n,int Gx_row,int Gz_row, GF2mat &Gx,GF2mat &Gz, GF2mat &Cx,GF2mat &Cz, int debug=0);





//files used in concatenated codes and product codes.

int reduce(GF2mat Gax, GF2mat Gaz, GF2mat Gbx, GF2mat Gbz,int ddax,int ddaz,int ddbx,int ddbz);


int concatenate(GF2mat Gax, GF2mat Gaz, GF2mat Gbx, GF2mat Gbz,int ddax,int ddaz,int ddbx,int ddbz);



// a version include both reduce and concatenation
// mode=1 for reduce/subsystem product
// mode=2 for concatenation
//only dz is checked cause dx is known to be tight
int product(GF2mat Gax, GF2mat Gaz, GF2mat Gbx, GF2mat Gbz,int ddax,int ddaz,int ddbx,int ddbz, int mode=1);




#endif
