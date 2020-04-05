//Weilei Zeng. Some small functions for use
#include <fstream>
#include<string>
#include<sstream>
#include "my_lib.h"
#include "lib.h"
#include <itpp/itbase.h>
using namespace itpp;
using namespace std;

bool is_quantum_code(GF2mat G_x, GF2mat G_z){
  //check if the code is quantum or not
  if (G_x.row_rank()+G_z.row_rank() < G_x.cols()){
    if ( (G_x*G_z.transpose()).is_zero()){
      return true;
    }
    cout<<"Not a quantum code:pertation relation is not satisfied."<<endl;
    return false;
  }else{
    cout<<"Not a quantum code: zero rank for codeword space."<<endl;
    return false;
  }    
}

GF2mat make_it_full_rank(GF2mat fat){
  //reduce a fat matrix with degenerate rows to a thin matrix with full rank; remove the dependent rows
  GF2mat thin, T,U;
  ivec P;
  int rank = fat.T_fact(T,U,P);
  if (rank = fat.rows()){ //if it is already full rank
    return fat;
  }
  thin = U.get_submatrix(0,0,rank-1,U.cols()-1);
  thin.permute_cols(P,true);
  return thin;
}

int GF2matPrint(GF2mat G,char * name){
  //print brief information of G
  cout<<name<<" ("<<G.rows()<<","<<G.cols()<<"), density = "
      <<G.density()<<endl;
  return 0;

}

GF2mat kron(GF2mat B, GF2mat A){//kronex tensor product of B x A
  //how could I define this in the wrong way?
  //although I define it in the wrong way, I dont think it will affect the result
  //but in any later calculation, use the right definition. May 26, 2018
  GF2mat G(1,1+A.cols()*B.cols());//an empty row; delete it finally
  GF2mat zero(A.rows(),A.cols());
  for (int i=0;i<B.rows();i++){
    GF2mat G_row(A.rows(),1);//an empty column, delete it finally
    for (int j=0;j<B.cols();j++){
      //cout<<"debug: j = "<<j<<endl;
      if (B.get(i,j)){//1,set A
	G_row = G_row.concatenate_horizontal(A);
      }else{//0, set zero
	G_row = G_row.concatenate_horizontal(zero);
      }
    }
    //    cout<<G<<endl<<G_row<<endl;
    G = G.concatenate_vertical(G_row);
  }
  
  //  cout<<G<<endl<<"debug"<<endl;
  G=G.get_submatrix(1,1,G.rows()-1,G.cols()-1);
  //cout<<G<<endl;
  return G;
}

string NumberToString(int pNumber) //not used anywhere. 
{
  ostringstream o;
  o<<pNumber;
  return o.str();
}

GF2mat append_vector(GF2mat G,bvec b){
  //append a row vector to the GF2mat G
  //used when saving error vectors for decoding
  G.set_size(G.rows()+1,G.cols(),true);
  G.set_row(G.rows()-1,b);
  return G;
}

GF2mat get_GF2mat(char * filename_prefix, char * filename_suffix){
  char filename[255];
  sprintf(filename,"%s%s",filename_prefix,filename_suffix);
  GF2mat E=MM_to_GF2mat(filename);
  return E;
}

GF2mat get_GF2mat(char * parent_folder, char * folder,  char *filename){
  char filename_full[255];
  sprintf(filename_full,"%s/%s/%s",parent_folder, folder, filename);
  GF2mat E=MM_to_GF2mat(filename_full);
  return E;
}

double get_error_density(GF2mat E){
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
