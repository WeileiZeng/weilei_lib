#ifndef MM_WRITE_H
#define MM_WRITE_H

//#include <stdio.h>
//#include <stdlib.h>
//#include <fstream> //ofstream ifstream
#include <itpp/itbase.h>
//using namespace itpp;
//using namespace std;

//#include <itpp/itbase.h>
//using namespace itpp;

//int GF2mat_to_MM(GF2mat G, char* file_name="mm_temp.dat");
int GF2mat_to_MM(itpp::GF2mat G, char* file_name, int debug=0);

//int mat_to_MM(mat G, char* file_name="mm_temp.dat");
int mat_to_MM(itpp::mat G, char* file_name);

//add string compatibility
int GF2mat_to_MM(itpp::mat G, std::string file_name);
int mat_to_MM(itpp::mat G, std::string file_name);

#endif

