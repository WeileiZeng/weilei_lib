//lib for general functions.


#ifndef LIB_H
#define LIB_H
#include <string>
#include <iostream>
#include<fstream>
#include <stdio.h>
#include <itpp/itbase.h>
using namespace itpp;
using namespace std;

bool is_quantum_code(GF2mat G_x, GF2mat G_z);
  //check if the code is quantum or not   
GF2mat make_it_full_rank(GF2mat fat);
  //reduce a fat matrix with degenerate rows to a thin matrix with full rank; remove the dependent rows 

int GF2matPrint(GF2mat G, char * name = (char *) " ");//print brief infomation of G
int matPrint(mat G, char * name = (char *) " ");//print brief infomation of G


GF2mat kron(GF2mat A, GF2mat B);//not sure how this kron is defined, maybe inversed

string NumberToString(int pNumber);//Convert number to string. //because the build in function is not supported in C++ 11


GF2mat append_vector(GF2mat G,bvec b);//add row vectors to a GF2mat

//getGF2mat from filename with given pattern
GF2mat get_GF2mat(char * filename_prefix, char * filename_suffix);
GF2mat get_GF2mat(char * parent_folder, char * folder, char * filename);

//fisrt row is zero
//remove first row and return density of the remained submatrix    
double get_error_density(GF2mat E);


//save mat into gnuplot data file, with a custom comment as header
int mat2gnudata(mat data, string filename, string header);


string color_text(string str);
//return red string. same as red_text

string red_text(string str);
//return red string

string blue_text(string str);
//return blue string


#endif
