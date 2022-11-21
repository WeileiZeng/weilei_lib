/**\file lib.cpp
 *\author Weilei Zeng
 *\brief implementation of lib.h
 */

//Weilei Zeng. Some small functions for use
//#include <fstream>
//#include<string>
//#include<iostream>
//#include<sstream>
//#include "my_lib.h"
//#include <stdio.h>
#include <itpp/itbase.h>
#include <chrono> //for time

#include "lib.h"
#include "mm_read.h"

/**
 * currently return red color as default; more color to be implemented
|     |    foreground| background|
|---|---|---|
|black        |30         |40|
|red          |31         |41|
|green        |32         |42|
|yellow       |33         |43|
|blue         |34         |44|
|magenta      |35         |45|
|cyan         |36         |46|
|white        |37         |47|
 */
std::string common::color_text(std::string str){
  //return text in red color
  return red_text(str);
}

std::string common::red_text(std::string str){
  //return text in red color
  return "\033[1;31m"+str+"\033[0m";
}

std::string common::blue_text(std::string str){
  //return text in blue color
  return "\033[1;34m"+str+"\033[0m";
}


/**
 * @param mode=1 secs (default)
 * @param mode=2 millisecs, \f$ 10 ^{-3} \f$
 * @param mode=3 10 nano secs = \f$10^{-8}\f$ sec
 */
int common::get_time(int mode){

  auto now = std::chrono::system_clock::now();
  //  auto now_ms = std::chrono::time_point_cast<std::chrono::milliseconds>(now);
  auto value = now.time_since_epoch();
  long duration = value.count();
  int t=0;
  const int DIGIT=1000000000;
  switch ( mode ){
  case 1: // seconds
    t = ( duration / 100000000 ) % DIGIT ;
    break;
  case 2: // milli seconds
    t = (duration / 100000) % DIGIT;
    break;
  case 3:
    t = duration % DIGIT;
    break;
  }
  return t;
}



bool common::is_quantum_code(itpp::GF2mat & G_x, itpp::GF2mat & G_z){
  //check if the code is quantum or not
  if (G_x.row_rank()+G_z.row_rank() < G_x.cols()){
    if ( (G_x*G_z.transpose()).is_zero()){
      return true;
    }
    std::cout<<"Not a quantum code:pertation relation is not satisfied."<<std::endl;
    return false;
  }else{
    std::cout<<"Not a quantum code: zero rank for codeword space."<<std::endl;
    return false;
  }    
}

bool common::is_quantum_code(itpp::GF2mat &Gx,itpp::GF2mat &Gz, itpp::GF2mat &Cx,itpp::GF2mat &Cz){
  if (!(Gx*Gz.transpose()).is_zero()){
    std::cout<<"(Gx*Gz.transpose()) is not zero"<<std::endl;
    //    throw "(Gx*Gz.transpose()) is not zero";
    return false;
  }
  if (!(Gx*Cz.transpose()).is_zero()){
    //    throw "(Gx*Cz.transpose()) is not zero";
    std::cout<<"(Gx*Cz.transpose()) is not zero"<<std::endl;return false;
  }
  if (!(Gz*Cx.transpose()).is_zero()){
    //    throw "(Gz*Cx.transpose()) is not zero";
    std::cout<<"(Gz*Cx.transpose()) is not zero"<<std::endl;return false;
  }
  int rank_of_Gx=Gx.row_rank();
  int rank_of_Gz=Gz.row_rank();
  int rank_of_Cx=Cx.row_rank();
  int rank_of_Cz=Cz.row_rank();
  int n=Gx.cols();
  if (rank_of_Gx+rank_of_Gz+rank_of_Cx != n){
    std::cout<<"(rank_of_Gx+rank_of_Gz+rank_of_Cx != n)"<<std::endl;
    return false;
  }
  if(rank_of_Cx != rank_of_Cz){
    std::cout<<"(rank_of_Cx != rank_of_Cz)"<<std::endl;return false;
  }
  //  std::cout<<"is_quantum_code(): It is a quantum code!"<<std::endl;
  return true;
}



itpp::GF2mat common::make_it_full_rank(itpp::GF2mat fat){
  //reduce a fat matrix with degenerate rows to a thin matrix with full rank; remove the dependent rows
  itpp::GF2mat thin, T,U;
  itpp::ivec P;
  int rank = fat.T_fact(T,U,P);
  if (rank == fat.rows()){ //if it is already full rank //this was only = but not ==, weilei changed on Dec 13, 2019, much later than the program was runned.
    return fat;
  }
  thin = U.get_submatrix(0,0,rank-1,U.cols()-1);
  thin.permute_cols(P,true);
  return thin;
}

/*
int itpp::GF2matPrint(itpp::GF2mat &G,char * name){
  //print brief information of G
  std::cout<<"GF2mat "<<name<<", size = ("<<G.rows()<<","<<G.cols()<<"), density = "
      <<G.density()<<std::endl;
  return 0;
}
*/
int common::GF2matPrint(itpp::GF2mat &G, std::string name){
  //print brief information of G
  //GF2matPrint(G, name.c_str());
    std::cout<<"GF2mat "<<name<<", size = ("<<G.rows()<<","<<G.cols()<<"), density = "
    <<G.density()<<std::endl;
  return 0;
}

/** should change char to string */
int common::matPrint(itpp::mat G,char * name){
  //print brief information of G
  std::cout<<"mat    "<<name<<", size = ("<<G.rows()<<","<<G.cols()<<")"<<std::endl;
  return 0;
}

/** 
 *\warning: The implementation is confusing, but should works as the conventional Kronecker product defined in the wikipedia page
 */
itpp::GF2mat common::kron(itpp::GF2mat A, itpp::GF2mat B){//Kroneker tensor product of A x B
  //how could I define this in the wrong way?
  //although I define it in the wrong way, I dont think it will affect the result
  //but in any later calculation, use the right definition. May 26, 2018

  //switch A B
  itpp::GF2mat  temp;
  temp = A;
  A=B;
  B=temp;
  
  itpp::GF2mat G(1,1+A.cols()*B.cols());//an empty row; delete it finally
  itpp::GF2mat zero(A.rows(),A.cols());
  for (int i=0;i<B.rows();i++){
    itpp::GF2mat G_row(A.rows(),1);//an empty column, delete it finally
    for (int j=0;j<B.cols();j++){
      //std::cout<<"debug: j = "<<j<<std::endl;
      if (B.get(i,j)){//1,set A
	G_row = G_row.concatenate_horizontal(A);
      }else{//0, set zero
	G_row = G_row.concatenate_horizontal(zero);
      }
    }
    //    std::cout<<G<<std::endl<<G_row<<std::endl;
    G = G.concatenate_vertical(G_row);
  }
  
  //  std::cout<<G<<std::endl<<"debug"<<std::endl;
  G=G.get_submatrix(1,1,G.rows()-1,G.cols()-1);
  //std::cout<<G<<std::endl;
  return G;
}

std::string common::NumberToString(int pNumber) //not used anywhere. 
{
  std::ostringstream o;
  o<<pNumber;
  return o.str();
}

itpp::GF2mat common::append_vector(itpp::GF2mat G,itpp::bvec b){
  //append a row vector to the GF2mat G
  //used when saving error vectors for decoding
  G.set_size(G.rows()+1,G.cols(),true);
  G.set_row(G.rows()-1,b);
  return G;
}

itpp::GF2mat common::get_GF2mat(char * filename_prefix, char * filename_suffix){
  char filename[255];
  sprintf(filename,"%s%s",filename_prefix,filename_suffix);
  itpp::GF2mat E=MM_to_GF2mat(filename);
  return E;
}

itpp::GF2mat common::get_GF2mat(char * parent_folder, char * folder,  char *filename){
  char filename_full[255];
  sprintf(filename_full,"%s/%s/%s",parent_folder, folder, filename);
  itpp::GF2mat E=MM_to_GF2mat(filename_full);
  return E;
}

double common::get_error_density(itpp::GF2mat E){
  //fisrt row is zero
  //remove first row and return density of the remained submatrix
  double w;
  if(E.rows()==1){
    w=0;
  }else{
    w=E.get_submatrix(1,0,E.rows()-1,E.cols()-1).density();
  }
  return w;
}

//save mat into gnuplot data file, with a custom comment as header
int common::mat2gnudata(itpp::mat data, std::string filename, std::string header){
  FILE *fout;//file to save the data
  //char filename_out[255];
  //sprintf(filename_out,"%s/gnuplot/rate_versus_p_size_%d_weight_%d.gnudat",error_folder,size,weight);
  fout =fopen(filename.c_str(),"w");
  fprintf(fout,"%s\n",header.c_str());
  int row=data.rows(), col=data.cols();
  for (int i =0;i<row;i++){
    for (int j=0;j<col;j++){
      fprintf(fout,"%f\t",data.get(i,j));
    }
    fprintf(fout,"\n");
  }
  fclose(fout);
  return 0;	      
}


/*
int common::weight(itpp::bvec &b){//we can use BPSK::count(zero_b(N),b) to count the weight
  int length=b.length();
  int weight=0;
  for (int i=0;i<length;i++){
    if(b.get(i)){
      weight++;
    }
  }
  return weight;
}
*/
