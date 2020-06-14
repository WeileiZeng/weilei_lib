//lib for general functions.


#ifndef LIB_H
#define LIB_H
#include <string>
#include <iostream>
#include<fstream>
#include <stdio.h>
#include <itpp/itbase.h>

//#include "weilei_lib.h"
//using namespace itpp;
//using namespace std;

bool is_quantum_code(itpp::GF2mat G_x, itpp::GF2mat G_z);
  //check if the code is quantum or not   
itpp::GF2mat make_it_full_rank(itpp::GF2mat fat);
  //reduce a fat matrix with degenerate rows to a thin matrix with full rank; remove the dependent rows 

//int itpp::GF2matPrint(itpp::GF2mat &G, char * name = (char *) " ");//print brief infomation of G
int GF2matPrint(itpp::GF2mat &G, std::string name);

int matPrint(itpp::mat G, char * name = (char *) " ");//print brief infomation of G


itpp::GF2mat kron(itpp::GF2mat A, itpp::GF2mat B);//not sure how this kron is defined, maybe inversed

std::string NumberToString(int pNumber);//Convert number to std::string. //because the build in function is not supported in C++ 11


itpp::GF2mat append_vector(itpp::GF2mat G, itpp::bvec b);//add row vectors to a itpp::GF2mat

//getitpp::GF2mat from filename with given pattern
itpp::GF2mat get_GF2mat(char * filename_prefix, char * filename_suffix);
itpp::GF2mat get_GF2mat(char * parent_folder, char * folder, char * filename);

//fisrt row is zero
//remove first row and return density of the remained submatrix    
double get_error_density(itpp::GF2mat E);


//save mat into gnuplot data file, with a custom comment as header
int mat2gnudata(itpp::mat data, std::string filename, std::string header);


std::string color_text(std::string str);
//return red std::string. same as red_text

std::string red_text(std::string str);
//return red std::string

std::string blue_text(std::string str);
//return blue std::string


int get_time(int mode = 1);

#endif
