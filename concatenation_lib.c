//outdated, not in use
//Weilei Zeng Nov 21, 2018
//to implement quantum concatenated codes. There are several ways of concatenation, see hypergraph_product_code.pdf


#include "dist.h"
#include <itpp/itbase.h>
#include <itpp/itcomm.h>
#include <stdio.h>
#include "my_lib.h"
#include "concatenation_lib.h"
using namespace itpp;
using namespace std;


int reduce(GF2mat Gax, GF2mat Gaz, GF2mat Gbx, GF2mat Gbz,int ddax,int ddaz,int ddbx,int ddbz){
  //construct code C and calculate the distance; Compare it with the input (estimated) value
  cout<<"estimate value of dax,daz,dbx,dbz = "<<ddax<<","<<ddaz<<","<<ddbx<<","<<ddbz<<","<<endl;
  int na=Gax.cols();//,nb=Gbx.cols();//,nc=na*nb;//size of the codes
  GF2mat Cax=getC(Gax,Gaz),Cbx=getC(Gbx,Gbz);//This line doesn't allow C to be empty
  GF2mat Caz=getC(Gax,Gaz,1),Cbz=getC(Gbx,Gbz,1);//This line doesn't allow C to be empty  
  GF2mat Gcz = kron(Gaz,Cbz).concatenate_vertical(kron(gf2dense_eye(na),Gbz));
  //  Gcz=make_it_full_rank(Gcz);//not sure if I need it here

  /*//check Cax
  if ( (Gaz*Cax.transpose()).is_zero() ){
      cout<<"Good Cax"<<endl;
      }*/
  GF2mat Gcx=kron(gf2dense_eye(na),Gbx).concatenate_vertical( kron(Gax,Cbx)   );
  //  Gcx=make_it_full_rank(Gcx);//not sure if I need it here
  /*  int daz = quantum_dist(Gax,Gaz,ddaz,1); //no need to check. Must euqal
  cout<<"daz="<<daz<<",ddaz="<<ddaz<<endl;
  int dbz = quantum_dist(Gbx,Gbz,ddbz,1);
  cout<<"dbz="<<dbz<<",ddbz="<<ddbz<<endl;*/
  /*  if (is_quantum_code(Gcx,Gcz) ){
    cout<<"C is a quantum Code."<<endl;
    }*/
  int daz=ddaz,dbz=ddbz;
  int dcz = quantum_dist(Gcx,Gcz,daz*dbz,1);//donot use estimated value ddaz and ddbz
  if (dcz == daz*dbz){
    cout<<"dcz = daz*dbz = "<<dcz<<endl;
  }else if(dcz == INF) {
    cout<<"dcz = "<<dcz<<", daz = "<<daz<<", dbz = "<<dbz<<endl;
  }else{
    cout<<"---------------------------------------------------------------------CASE: daz*dbz="<<daz*dbz<<", dcz="<<dcz<<endl;
  }
  return 0;
  
}

int concatenate(GF2mat Gax, GF2mat Gaz, GF2mat Gbx, GF2mat Gbz,int ddax,int ddaz,int ddbx,int ddbz){
  //construct code C and calculate the distance; Compare it with the input (estimated) value
  cout<<"estimate value of dax,daz,dbx,dbz = "<<ddax<<","<<ddaz<<","<<ddbx<<","<<ddbz<<","<<endl;
  int na=Gax.cols(),nb=Gbx.cols();//,nc=na*nb;//size of the codes
  GF2mat Gcz = kron(Gaz,gf2dense_eye(nb)).concatenate_vertical(kron(gf2dense_eye(na),Gbz));
  //  Gcz=make_it_full_rank(Gcz);//not sure if I need it here
  GF2mat Cax=getC(Gax,Gaz),Cbx=getC(Gbx,Gbz);//This line doesn't allow C to be empty
  /*//check Cax
  if ( (Gaz*Cax.transpose()).is_zero() ){
      cout<<"Good Cax"<<endl;
      }*/
  GF2mat Gcx=kron(Gax,Gbx).concatenate_vertical(
       	         kron(Cax,Gbx).concatenate_vertical(
          	     kron(Gax,Cbx)
								  )
						);
  //  Gcx=make_it_full_rank(Gcx);//not sure if I need it here
  /*  int daz = quantum_dist(Gax,Gaz,ddaz,1);
  cout<<"daz="<<daz<<",ddaz="<<ddaz<<endl;
  int dbz = quantum_dist(Gbx,Gbz,ddbz,1);
  cout<<"dbz="<<dbz<<",ddbz="<<ddbz<<endl;*/
  /*  if (is_quantum_code(Gcx,Gcz) ){
    cout<<"C is a quantum Code."<<endl;
    }*/
  int daz=ddaz,dbz=ddbz;
  int dcz = quantum_dist(Gcx,Gcz,daz*dbz,1);//donot use estimated value ddaz and ddbz
  if (dcz == daz*dbz){
    cout<<"dcz = daz*dbz = "<<dcz<<endl;
  }else if(dcz == INF) {
    cout<<"dcz = "<<dcz<<", daz = "<<daz<<", dbz = "<<dbz<<endl;
  }else{
    cout<<"---------------------------------------------------------------------CASE: daz*dbz="<<daz*dbz<<", dcz="<<dcz<<endl;
  }
  return 0;
  
}
