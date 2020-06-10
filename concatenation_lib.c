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



bool is_quantum_code(GF2mat &Gx,GF2mat &Gz, GF2mat &Cx,GF2mat &Cz){
  if (!(Gx*Gz.transpose()).is_zero()){
    cout<<"(Gx*Gz.transpose()) is not zero"<<endl;
    //    throw "(Gx*Gz.transpose()) is not zero";
    return false;
  }
  if (!(Gx*Cz.transpose()).is_zero()){
    //    throw "(Gx*Cz.transpose()) is not zero";
    cout<<"(Gx*Cz.transpose()) is not zero"<<endl;return false;
  }
  if (!(Gz*Cx.transpose()).is_zero()){
    //    throw "(Gz*Cx.transpose()) is not zero";
    cout<<"(Gz*Cx.transpose()) is not zero"<<endl;return false;
  }
  int rank_of_Gx=Gx.row_rank();
  int rank_of_Gz=Gz.row_rank();
  int rank_of_Cx=Cx.row_rank();
  int rank_of_Cz=Cz.row_rank();
  int n=Gx.cols();
  if (rank_of_Gx+rank_of_Gz+rank_of_Cx != n){
    cout<<"(rank_of_Gx+rank_of_Gz+rank_of_Cx != n)"<<endl;
    return false;
  }
  if(rank_of_Cx != rank_of_Cz){
    cout<<"(rank_of_Cx != rank_of_Cz)"<<endl;return false;
  }
  //  cout<<"is_quantum_code(): It is a quantum code!"<<endl;
  return true;
}

/*int getRandomQuantumCode(GF2mat &Gx,GF2mat &Gz, GF2mat &Cx,GF2mat &Cz){
  int n=21;//sample input
  int Gx_row=8;
  int Gz_row=8;
  getRandomQuantumCode(n,Gx_row,Gz_row,Gx,Gz,Cx,Cz);
  return 0;
  }*/

GF2mat remove_col(GF2mat G, int col){
  int n = G.cols();
  if ( col ==0 ) return G.get_submatrix(0,1,G.rows()-1,G.cols()-1);
  if ( col == n-1 ) return G.get_submatrix(0,0,G.rows()-1,G.cols()-2);
  return G.get_submatrix(0,0,G.rows()-1,col-1).concatenate_horizontal(
								      G.get_submatrix(0,col+1,G.rows()-1,G.cols()-1)
								      );  
}

void remove_singleton(GF2mat &Gx, GF2mat &Gz){
  //remove zero columns in Gx and Gz
  //not in use, just discard code with distance 1, easier solution
  int n = Gx.cols();
  bvec to_remove(n);//1 fro remove, 0 remain  
  for ( int i=0;i<n;i++){
    if (GF2mat(Gx.get_col(i)).is_zero()) to_remove.set(i,1);
    if (GF2mat(Gz.get_col(i)).is_zero()) to_remove.set(i,1);
  }
  for ( int i=0;i<n;i++){
    if ( to_remove(n-i-1) ){
	Gx=remove_col(Gx,n-i-1);
	Gz=remove_col(Gz,n-i-1);
      }
  }
  Gx=make_it_full_rank(Gx);
  Gz=make_it_full_rank(Gz);
  cout<<"singleton removed"<<endl;
  return;
}

int getRandomQuantumCode(int n,int Gx_row,int Gz_row, GF2mat &Gx,GF2mat &Gz, GF2mat &Cx,GF2mat &Cz){

  Gx = GF2mat(Gx_row,n);
  Gz = GF2mat(Gz_row,n);
  for ( int i =0;i<Gx_row;i++){//random G_x
    Gx.set_row(i,randb(n));//equally 0 and 1s
  }
  //Gx might not be full rank at this point
  
  GF2mat T,U; ivec P;
  int rank_of_Gx = Gx.transpose().T_fact(T,U,P);
  GF2mat Q=T.get_submatrix(rank_of_Gx,0,n-1,n-1);


  //  if ( debug ) 
  //  cout<<"Gx 1st row:"<<Gx.get_row(0)<<endl; // for debug the random seed

  //  if ( rank_of_Gx < Gx.rows() ) cout<<"getRandomQuantumCode: Gx not full rank"<<endl;
  //else cout<<"getRandomQuantumCode: Gx is  full rank"<<endl;
  
  //  Q.permute_cols(P,true); no need for T, only need for U which is not used here
  //  GF2matPrint(Q,"Q");
  GF2mat alpha(Gz_row,Q.rows()); //a random binary matrix to select G_z
  for ( int i=0;i<Gz_row;i++){
    alpha.set_row(i,randb(Q.rows()));
  }
  Gz=alpha*Q;
  //  Gz=Q.get_submatrix(0,0,Gz_row-1,n-1);
  //Cz=Q.get_submatrix(Gz_row,0,Q.rows()-1,n-1);
  //  GF2matPrint(Gz,"Gz");
  //  GF2matPrint(Cz,"Cz");
  // the following 2 are bad trials, which make distance 1 always
  //Gx = nullSpace(Gz).get_submatrix(0,0,Gx_row-1,n-1);
  //Gz = nullSpace(Gx).get_submatrix(0,0,Gz_row-1,n-1);
  //  remove_singleton(Gx,Gz);
  Cx=getC(Gx,Gz);
  Cz=getC(Gx,Gz,1);
  //  if (! is_quantum_code(Gx,Gz,Cx,Cz)) throw "invalid code";

  return 0;
}

int getGoodQuantumCode(int n,int Gx_row,int Gz_row, GF2mat &Gx,GF2mat &Gz, GF2mat &Cx,GF2mat &Cz, int debug){
  // return best codes among multip trial
  //repeat multiple times to get the best distance
  GF2mat Gx_temp, Gz_temp,Cx_temp,Cz_temp;
  int search_trial=1000;
  int flag_find_good_code=0;
  for ( int i =0; i<search_trial; i++){

    getRandomQuantumCode( n, Gx_row,Gz_row, Gx_temp, Gz_temp,Cx_temp,Cz_temp);
    //check distance and update if get larger distance
    int dx = quantum_dist_v2(Gx_temp,Gz_temp);
    if ( dx >1 ){
      int dz = quantum_dist_v2(Gx_temp,Gz_temp,1);
      if (dz >1 ){
	flag_find_good_code=1;
	//	Gx = Gx_temp; Gz = Gz_temp; Cx = Cx_temp; Cz = Cz_temp;
	if (debug) cout<<blue_text("get good code when i =")<<i<<endl;
	break;
	//	return 0;
      }
    }
  }
  Gx = Gx_temp; Gz = Gz_temp; Cx = Cx_temp; Cz = Cz_temp;

  //  if ( flag_find_good_code){
  if ( debug ) cout<<"Gx 1st row:"<<Gx.get_row(0)<<endl; // for debug the random seed
  if ( Gx.row_rank() < Gx.rows() ) {
    if (debug) cout<<"getGoodQuantumCode: Gx not full rank. now make it full rank"<<endl;
    Gx = make_it_full_rank(Gx);
  }

  if ( debug) if ( ! flag_find_good_code ) cout<<color_text("didn't find good code after ")<<search_trial<<" trials"<<endl;
  return 0;
}

void set_submatrix(GF2mat & G, GF2mat sub, int row, int col){
  //put sub into G, start from (row,col)
  for ( int i =0 ; i < sub.rows(); i ++)
    for ( int j = 0; j< sub.cols(); j++ )
      G.set(i+row, j+col, sub.get(i,j));
  return;
}

int check(){
  return 0;
}

// generate all code with size na systematically
int generate_code(GF2mat & Gax, GF2mat & Gaz, int na, int Gax_row, int id_Gax, int id_Gaz){
  Gax = GF2mat(Gax_row,na);
  // identity matrix in the left part to make it reduce row echelon form.
  for ( int i =0;i<Gax_row; i++) Gax.set(i,i,1) ;
  //id in (0,2^( (Gax_row * (na-Gax_row) ))
  // check id
  const id_Gax_MAX = (int) pow(2,  Gax_row * (na-Gax_row) ) -2 ;
  if ( id_Gax <1 || id_Gax > id_Gax_MAX -1 ) {
    cout<<"illegal id_Gax"<<endl;
    throw 2;
  }
  GF2mat alpha_Gax = GF2mat( dec2bin(Gax_row*(na-Gax_row), i), false);//false for row vector
  for ( int i = 0 ; i < Gax_row; i++)
    set_submatrix(Gax,alpha_Gax.get_submatrix(0, i*(na-Gax_row), 0, (i+1)*(na-Gax_row)-1), i, na-Gax_row);
  GF2mat H = nullSpace(alpha_Gax);
  //check id_Gaz
  const id_Gaz_MAX = (int) pow(2, na - Gax_row) - 2;
  if ( id_Gaz < 1 || id_Gaz > id_Gaz_MAX -1 ){
    cout<<"illegal id_Gaz"<<endl;
    throw 2;
  }
  GF2mat alpha_Gaz = GF2mat(dec2bin(id_Gaz), false);
  Gaz = alpha_Gaz*H;
  
  return 0;
}





/*

int concatenate(GF2mat Gax, GF2mat Gaz, GF2mat Gbx, GF2mat Gbz,int ddax,int ddaz,int ddbx,int ddbz){
  //construct code C and calculate the distance; Compare it with the input (estimated) value
  cout<<"estimate value of dax,daz,dbx,dbz = "<<ddax<<","<<ddaz<<","<<ddbx<<","<<ddbz<<","<<endl;
  int na=Gax.cols();//,nb=Gbx.cols();//,nc=na*nb;//size of the codes
  GF2mat Cax=getC(Gax,Gaz),Cbx=getC(Gbx,Gbz);//This line doesn't allow C to be empty
  GF2mat Caz=getC(Gax,Gaz,1),Cbz=getC(Gbx,Gbz,1);//This line doesn't allow C to be empty  
  GF2mat Gcz = kron(Gaz,Cbz).concatenate_vertical(kron(gf2dense_eye(na),Gbz));
  //  Gcz=make_it_full_rank(Gcz);//not sure if I need it here


  GF2mat Gcx=kron(gf2dense_eye(na),Gbx).concatenate_vertical( kron(Gax,Cbx)   );
  //  Gcx=make_it_full_rank(Gcx);//not sure if I need it here

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

int reduce(GF2mat Gax, GF2mat Gaz, GF2mat Gbx, GF2mat Gbz,int ddax,int ddaz,int ddbx,int ddbz){
  //construct code C and calculate the distance; Compare it with the input (estimated) value
  cout<<"estimate value of dax,daz,dbx,dbz = "<<ddax<<","<<ddaz<<","<<ddbx<<","<<ddbz<<","<<endl;
  int na=Gax.cols(),nb=Gbx.cols();//,nc=na*nb;//size of the codes
  GF2mat Gcz = kron(Gaz,gf2dense_eye(nb)).concatenate_vertical(kron(gf2dense_eye(na),Gbz));
  //  Gcz=make_it_full_rank(Gcz);//not sure if I need it here
  GF2mat Cax=getC(Gax,Gaz),Cbx=getC(Gbx,Gbz);//This line doesn't allow C to be empty

  GF2mat Gcx=kron(Gax,Gbx).concatenate_vertical(
       	         kron(Cax,Gbx).concatenate_vertical(
          	     kron(Gax,Cbx)
								  )
						);
  //  Gcx=make_it_full_rank(Gcx);//not sure if I need it here

  int daz=ddaz,dbz=ddbz;
  int dcz = quantum_dist(Gcx,Gcz,daz*dbz,1);//donot use estimated value ddaz and ddbz
  if (dcz == daz*dbz){
    cout<<"dcz = daz*dbz = "<<dcz<<endl;
  }else if(dcz == INF) {
    cout<<"dcz = "<<dcz<<", daz = "<<daz<<", dbz = "<<dbz<<endl;
  }else{
    cout<<"-------------------------------------------------------------counter example CASE: daz*dbz="<<daz*dbz<<", dcz="<<dcz<<endl;
  }
  return 0;
  
}
*/


// a version include both reduce and concatenation
// mode=1 for reduce/subsystem product
// mode=2 for concatenation
//only dz is checked cause dx is known to be tight
int product(GF2mat Gax, GF2mat Gaz, GF2mat Gbx, GF2mat Gbz,int ddax,int ddaz,int ddbx,int ddbz, int debug, int mode){
  //construct code C and calculate the distance; Compare it with the input (estimated) value
  int na=Gax.cols(),nb=Gbx.cols();//,nc=na*nb;//size of the codes

  GF2mat Cax=getC(Gax,Gaz),Cbx=getC(Gbx,Gbz);//This line doesn't allow C to be empty
  GF2mat Caz=getC(Gax,Gaz,1),Cbz=getC(Gbx,Gbz,1);//This line doesn't allow C to be empty  
  
  //  Gcz=make_it_full_rank(Gcz);//not needed for calculating distance
  GF2mat Gcx,Gcz;

  switch ( mode ){
  case 0://reduce/subsystem product, x distance
    
    Gcx = kron(Gax,gf2dense_eye(nb)).concatenate_vertical(kron(gf2dense_eye(na),Gbx));
    Gcz=kron(Gaz,Gbz).concatenate_vertical(
					   kron(Caz,Gbz)
					   .concatenate_vertical(kron(Gaz,Cbz))
					   );
    break;
  case 1://reduce/subsystem product, z distance
    Gcz = kron(Gaz,gf2dense_eye(nb)).concatenate_vertical(kron(gf2dense_eye(na),Gbz));
    Gcx=kron(Gax,Gbx)
      .concatenate_vertical(
			    kron(Cax,Gbx)
			    .concatenate_vertical(kron(Gax,Cbx))
			    );
    break;
  case 2://concatenation
    Gcz = kron(Gaz,Cbz).concatenate_vertical(kron(gf2dense_eye(na),Gbz));
    Gcx=kron(gf2dense_eye(na),Gbx).concatenate_vertical( kron(Gax,Cbx)   );
    break;
  case 3:
    // chain complex to two CSS codes.
  case 4:
    {
    // chain complex to two CSS codes.
    Gcx=kron(Gaz.transpose(), gf2dense_eye(Gbx.rows()))
      .concatenate_horizontal(kron(gf2dense_eye(Gax.cols()),Gbx))
      .concatenate_horizontal(kron(GF2mat(Gax.cols(),Gax.rows()),GF2mat(Gbx.rows(),Gbz.rows())));
    Gcx = Gcx
      .concatenate_vertical(
			    kron(GF2mat(Gax.rows(),Gaz.rows()),GF2mat(Gbx.cols(),Gbx.rows()))
			    .concatenate_horizontal(kron(Gax, gf2dense_eye(Gbx.cols())))
			    .concatenate_horizontal(kron(gf2dense_eye(Gax.rows()),Gbz.transpose()))
			    );
    Gcz=kron(gf2dense_eye(Gaz.rows()), Gbx.transpose())
      .concatenate_horizontal(kron(Gaz,gf2dense_eye(Gbx.cols())))
      .concatenate_horizontal(kron(GF2mat(Gaz.rows(),Gax.rows()),GF2mat(Gbz.cols(),Gbz.rows())));
    Gcz=Gcz
      .concatenate_vertical(
			    kron(GF2mat( Gaz.cols(),Gaz.rows() ), GF2mat( Gbz.rows(),Gbx.rows() ))
			    .concatenate_horizontal(kron(gf2dense_eye(Gaz.cols()),Gbz))
			    .concatenate_horizontal(kron(Gax.transpose(),gf2dense_eye(Gbz.rows())))
			    );
    
    }
    break;

  }  


  int flag_dist_flip=1;
  switch ( mode ){
  case 0:
  case 4:
    //x distance
    flag_dist_flip=0;
    break;
  case 1:
  case 2:
  case 3:
    // z distance
    flag_dist_flip=1;
    break;
  }
  if ( debug ){
    switch ( mode ){
    case 3:
      cout<<"mode (3)"<<endl;
      break;
    case 4:
      cout<<"mode (4)"<<endl;
      break;
    }
  }
  if (debug){GF2matPrint(Gcx,"Gcx"); GF2matPrint(Gcz,"Gcz");}

  if ( ! (Gcx*Gcz.transpose()).is_zero() ){
    cout<<"concatenation_lib: not a quantum code "<<endl;
    throw "not a quantum code";
  }else{
    if ( debug) cout<<"mode ("<<mode<<") is quantum code"<<endl;
  }

  //  Gcx=make_it_full_rank(Gcx);//not sure if I need it here
  /*  int daz = quantum_dist(Gax,Gaz,ddaz,1);
  cout<<"daz="<<daz<<",ddaz="<<ddaz<<endl;
  int dbz = quantum_dist(Gbx,Gbz,ddbz,1);
  cout<<"dbz="<<dbz<<",ddbz="<<ddbz<<endl;*/
  /*  if (is_quantum_code(Gcx,Gcz) ){
    cout<<"C is a quantum Code."<<endl;
    }*/
  //  if ( debug ) cout<<"Gcx"<<Gcx<<"Gcz"<<Gcz<<endl;

  switch ( flag_dist_flip ){
  case 0:
    //x distance
    {
      //  if ( debug ) cout<<"Gcx"<<Gcx<<"Gcz"<<Gcz<<endl;
      int dax=ddax,dbx=ddbx;
      int dcx = quantum_dist(Gcx,Gcz,dax*dbx,debug,0);//donot use estimated value ddaz and ddbz
      if (debug) cout<<"dax,daz,dbx,dbz = "<<ddax<<","<<ddaz<<","<<ddbx<<","<<ddbz<<","<<endl;    
      if (dcx == dax*dbx){
	if (debug) cout<<"dcx = dax*dbx = "<<dcx<<endl;
	return 0;
      }else if(dcx == INF) {
	if (debug) cout<<"dcx = "<<dcx<<", dax = "<<dax<<", dbx = "<<dbx<<endl;
	return 1;
      }else{
	if (dcx > dax*dbx) cout<<"PSEUDO ";
	cout<<red_text("CASE:")<<" mode ("<<mode<<") dax*dbx="<<dax*dbx<<", dcx="<<dcx;
	cout<<". dax,daz,dbx,dbz = "<<ddax<<","<<ddaz<<","<<ddbx<<","<<ddbz<<";";    
	cout<<"na,nb,nc,"<<Gax.cols()<<","<<Gbx.cols()<<","<<Gcx.cols()<<";";
	cout<<"ka,kb="<<Cax.rows()<<","<<Cbx.rows()<<";";
	return 2;
      }
    }

    break;
  case 1:
    //z distance
    {
      int daz=ddaz,dbz=ddbz;
      int dcz = quantum_dist(Gcx,Gcz,daz*dbz,debug,1);//donot use estimated value ddaz and ddbz
      if (debug) cout<<"dax,daz,dbx,dbz = "<<ddax<<","<<ddaz<<","<<ddbx<<","<<ddbz<<","<<endl;    
      if (dcz == daz*dbz){
	if (debug) cout<<"dcz = daz*dbz = "<<dcz<<endl;
	return 0;
      }else if(dcz == INF) {
	if (debug) cout<<"dcz = "<<dcz<<", daz = "<<daz<<", dbz = "<<dbz<<endl;
	return 1;
      }else{
	if (dcz > daz*dbz) cout<<"PSEUDO ";
	cout<<red_text("CASE:")<<" mode ("<<mode<<") daz*dbz="<<daz*dbz<<", dcz="<<dcz;
	cout<<". dax,daz,dbx,dbz = "<<ddax<<","<<ddaz<<","<<ddbx<<","<<ddbz<<";";
	cout<<"na,nb,nc="<<Gax.cols()<<","<<Gbx.cols()<<","<<Gcx.cols()<<";";    
	cout<<"ka,kb="<<Cax.rows()<<","<<Cbx.rows()<<";";
	return 2;
      }
    }
  }
  return 0;
  
}
