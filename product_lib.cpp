//outdated, not in use
//Weilei Zeng Nov 21, 2018
//to implement quantum concatenated codes. There are several ways of concatenation, see hypergraph_product_code.pdf



#include "dist.h"
#include <itpp/itbase.h>
#include <itpp/itcomm.h>
#include <stdio.h>
//#include "my_lib.h"
//#include "concatenation_lib.h"
#include "product_lib.h"
#include<typeinfo> //for typeid(a).name()
//using namespace common;



/** Print information for CSS code*/
template <class CodeType>
std::ostream& print_code(std::ostream& os, const CodeType& code){
  os<<code.title<<"("<< code.type<<"): "
    <<"[n,k,dx,dz]=["<<code.n<<","<<code.k<<","<<code.dx<<","<<code.dz<<"]";
  return os;
}

std::ostream& operator<<(std::ostream& os, const ClassicalCode& code){
  //  os<<"---"<< code.type<<":"<<code.title
  os<<code.title<<"("<< code.type<<"): "
    <<" [n,k,d]=["<<code.n<<","<<code.k<<","<<code.d<<"]";
  return os;
  //  return print_code(os, code);
}
std::ostream& operator<<(std::ostream& os, const CSSCode& code){
  return print_code(os, code);
}
std::ostream& operator<<(std::ostream& os, const ProductCSSCode& code){
  return print_code(os, code);
}
std::ostream& operator<<(std::ostream& os, const SubsystemProductCSSCode& code){
  return print_code(os, code);
}
std::ostream& operator<<(std::ostream& os, const ConcatenatedProductCSSCode& code){
  return print_code(os, code);
}




  //constructor
ClassicalCode::ClassicalCode(){
}
ClassicalCode::ClassicalCode(itpp::GF2mat G, itpp::GF2mat H){
  G=G;H=H;
  return;
}

//distance estimator
int ClassicalCode::dist(){
  return rand_dist();
}
int ClassicalCode::min_weight_dist(){
  return common::min_wt_decoding(G);
}
int ClassicalCode::rand_dist(){
  return common::rand_dist(G);
}



//function
void ClassicalCode::info(){
  std::cout<<"Classical code: n = "<<n<<std::endl
	   <<"codeword generating matrix G"<<G
	   <<"parity check matrix H"<<H<<std::endl;
  return;
}
ClassicalCode ClassicalCode::dual(){
  ClassicalCode dual_code(H,G);
  return dual_code;
}

void ClassicalCode::full_rank(){
  G = common::make_it_full_rank(G);
  H = common::make_it_full_rank(H);
  n = G.cols();
  k = G.rows();
  if (k + H.rows() != n) {
    std::cout<<"ClassicalCode: This code is not valid"<<std::endl;
    throw 2;
  }
  return;
}

//generate sample code
void ClassicalCode::get_repetition_code(int L){
  n=L;
  H =  common::get_check_rept(n);
  G = itpp::GF2mat(itpp::ones_b(n), false);
  return;
}

void ClassicalCode::get_743_code(int L){
  n = L;
  H =  common::get_check_code743(n);
  G =  common::get_check_code734(n);
  //  std::cout<<"check \n"<<G<<H;
  return;
}

void ClassicalCode::get_734_code(int L){
  n=L;
  H =  common::get_check_code734(n);
  G =  common::get_check_code743(n);
  return;
}



CSSCode::CSSCode(){}
CSSCode::CSSCode(int na, int Gax_row, int id_Gax, int Gaz_row, int id_Gaz){
    n=na;Gx_row=Gax_row; id_Gx=id_Gax;
    Gz_row=Gaz_row;id_Gz=id_Gaz;
    is_defined=1;
}
int CSSCode::generate_by_id(int debug){

  int temp = generate_code(Gx, Gz, n, Gx_row, id_Gx, Gz_row, id_Gz, debug);
  if ( temp == 2){
    if ( true ) std::cout<<"Duplicate code for this ID, code not generated"<<std::endl;
  }
  return temp;
}

/** the G matrices may not be full rank*/
int CSSCode::getRandomCode(){
  return getRandomQuantumCode(n, Gx_row, Gz_row, Gx, Gz, Cx, Cz);
}

int CSSCode::getGoodCode(int debug){
  return getGoodQuantumCode(n, Gx_row, Gz_row, Gx, Gz, Cx, Cz, debug);
}


int CSSCode::set_up_CxCz(){  
  Cx=common::getC(Gx,Gz);
  Cz=common::getC(Gz,Gx);  
  is_C_defined=1;
  return 0;
}

bool CSSCode::is_valid(){
  if ( is_C_defined ) {
    return common::is_quantum_code(Gx, Gz, Cx, Cz);
  }
  return common::is_quantum_code(Gx, Gz);
}

void CSSCode::full_rank(){
  Gx = common::make_it_full_rank(Gx);
  Gz = common::make_it_full_rank(Gz);
  Cx = common::make_it_full_rank(Cx);
  Cz = common::make_it_full_rank(Cz);
  n = Gx.cols();
  k = Cx.rows();
  if (Gz.cols()==n && Cx.cols()==n && Cz.cols()==n \
      && Cz.rows()==k \
      && Gx.rows()+Gz.rows()+k ==n){
  }else{
    std::cout<<"CSSCode: This code is not valid"<<std::endl;
    throw 2;
  }
  return;
}

/** operator<< is better than this function*/
void CSSCode::info(){
  //  std::cout<<"info()"<<this<<std::endl;
  std::cout<<"info():"<<"[n,k,dx,dz]=["<<n<<","<<k<<","<<dx<<","<<dz<<"]"<<std::endl;
  return;
}

void CSSCode::dist(){
  dx = rand_dist_x();
  dz = rand_dist_z();
  return;
}

int CSSCode::min_weight_dist_x(){
  return  common::min_wt_decoding(Cx, Gx);
}
int CSSCode::min_weight_dist_z(){
  return  common::min_wt_decoding(Cz, Gz);
}

int CSSCode::rand_dist_x(){
  return common::quantum_dist_v2(Gx, Gz);
}
int CSSCode::rand_dist_z(){
  return common::quantum_dist_v2(Gx, Gz, 1);
}

void decode(itpp::bvec e_in, itpp::bvec e_out){


}

void CSSCode::get_713_code(){
  Gx=common::get_check_code743(n);
  Gz=common::get_check_code743(n);
  return;

}

int CSSCode::save(std::string filename_prefix){
  const char * title = filename_prefix.c_str();
  char filename_Gx[256];char filename_Gz[256];
  sprintf(filename_Gx,"%sGx.mm",title);  
  sprintf(filename_Gz,"%sGz.mm",title); 
  GF2mat_to_MM(Gx,filename_Gx);
  GF2mat_to_MM(Gz,filename_Gz);
  return 0;
}

int CSSCode::load(std::string filename_prefix){
  Gx=MM_to_GF2mat(filename_prefix+"Gx.mm");
  Gz=MM_to_GF2mat(filename_prefix+"Gz.mm");
  n=Gx.cols();  
  return 0;
}



ProductCSSCode::ProductCSSCode(CSSCode codeA_temp, CSSCode codeB_temp){
  //  SubsystemProductCSSCode(CSSCode codeA_temp, CSSCode codeB_temp){
    codeA=codeA_temp;
    codeB=codeB_temp;
    //    std::cout<<" get codeA with codeA.n = "<< codeA_temp.n<<std::endl;
    if ( codeA.is_defined && codeB.is_defined ){
      //      std::cout<<"both code A and code B are defined"<<std::endl;
      is_defined=1;
    }
  
}

void SubsystemProductCSSCode::product(){
  Gx = common::kron(codeA.Gx,itpp::gf2dense_eye(codeB.n)).concatenate_vertical(common::kron(itpp::gf2dense_eye(codeA.n),codeB.Gx));
  Gz = common::kron(codeA.Gz,itpp::gf2dense_eye(codeB.n)).concatenate_vertical(common::kron(itpp::gf2dense_eye(codeA.n),codeB.Gz));
  
  /*  std::cout<<"debug 0:"<<std::endl;
  std::cout<<codeA.Gz<<codeB.Gz<<std::endl;
  Hx = common::kron(codeA.Gz,codeB.Gz);
  std::cout<<"debug 0.1:"<<std::endl; */
  Hz=common::kron(codeA.Gz,codeB.Gz)
    .concatenate_vertical(common::kron(codeA.Cz,codeB.Gz))
    .concatenate_vertical(common::kron(codeA.Gz,codeB.Cz));
  //  std::cout<<"debug 1:"<<std::endl;
  Hx=common::kron(codeA.Gx,codeB.Gx)
    .concatenate_vertical(common::kron(codeA.Cx,codeB.Gx))
    .concatenate_vertical(common::kron(codeA.Gx,codeB.Cx));

 
  return;
}

/*SubsystemProductCSSCode::SubsystemProductCSSCode(CSSCode codeA_temp, CSSCode codeB_temp):ProductCSSCode::ProductCSSCode(CSSCode codeA_temp, CSSCode codeB_temp){
  }*/



itpp::GF2mat remove_col(itpp::GF2mat G, int col){
  int n = G.cols();
  if ( col ==0 ) return G.get_submatrix(0,1,G.rows()-1,G.cols()-1);
  if ( col == n-1 ) return G.get_submatrix(0,0,G.rows()-1,G.cols()-2);
  return G.get_submatrix(0,0,G.rows()-1,col-1).concatenate_horizontal(
								      G.get_submatrix(0,col+1,G.rows()-1,G.cols()-1)
								      );  
}

void remove_singleton(itpp::GF2mat &Gx, itpp::GF2mat &Gz){
  //remove zero columns in Gx and Gz
  //not in use, just discard code with distance 1, easier solution
  int n = Gx.cols();
  itpp::bvec to_remove(n);//1 fro remove, 0 remain  
  for ( int i=0;i<n;i++){
    if (itpp::GF2mat(Gx.get_col(i)).is_zero()) to_remove.set(i,1);
    if (itpp::GF2mat(Gz.get_col(i)).is_zero()) to_remove.set(i,1);
  }
  for ( int i=0;i<n;i++){
    if ( to_remove(n-i-1) ){
	Gx=remove_col(Gx,n-i-1);
	Gz=remove_col(Gz,n-i-1);
      }
  }
  Gx=common::make_it_full_rank(Gx);
  Gz=common::make_it_full_rank(Gz);
  std::cout<<"singleton removed"<<std::endl;
  return;
}

int getRandomQuantumCode(int n,int Gx_row,int Gz_row, itpp::GF2mat &Gx,itpp::GF2mat &Gz, itpp::GF2mat &Cx,itpp::GF2mat &Cz){

  Gx = itpp::GF2mat(Gx_row,n);
  Gz = itpp::GF2mat(Gz_row,n);
  for ( int i =0;i<Gx_row;i++){//random G_x
    Gx.set_row(i,itpp::randb(n));//equally 0 and 1s
  }
  //Gx might not be full rank at this point
  
  itpp::GF2mat T,U; itpp::ivec P;
  int rank_of_Gx = Gx.transpose().T_fact(T,U,P);
  itpp::GF2mat Q=T.get_submatrix(rank_of_Gx,0,n-1,n-1);


  itpp::GF2mat alpha(Gz_row,Q.rows()); //a random binary matrix to select G_z
  for ( int i=0;i<Gz_row;i++){
    alpha.set_row(i,itpp::randb(Q.rows()));
  }
  Gz=alpha*Q;
  Cx=common::getC(Gx,Gz);
  Cz=common::getC(Gx,Gz,1);
  //  if (! is_quantum_code(Gx,Gz,Cx,Cz)) throw "invalid code";
  return 0;
}

int getGoodQuantumCode(int n,int Gx_row,int Gz_row, itpp::GF2mat &Gx,itpp::GF2mat &Gz, itpp::GF2mat &Cx,itpp::GF2mat &Cz, int debug){
  // return best codes among multip trial
  //repeat multiple times to get the best distance
  itpp::GF2mat Gx_temp, Gz_temp,Cx_temp,Cz_temp;
  int search_trial=1000;
  int flag_find_good_code=0;
  for ( int i =0; i<search_trial; i++){

    getRandomQuantumCode( n, Gx_row,Gz_row, Gx_temp, Gz_temp,Cx_temp,Cz_temp);
    //check distance and update if get larger distance
    int dx = common::quantum_dist_v2(Gx_temp,Gz_temp);
    if ( dx >1 ){
      int dz = common::quantum_dist_v2(Gx_temp,Gz_temp,1);
      if (dz >1 ){
	flag_find_good_code=1;
	//	Gx = Gx_temp; Gz = Gz_temp; Cx = Cx_temp; Cz = Cz_temp;
	if (debug) std::cout<<common::blue_text("get good code when i =")<<i<<std::endl;
	break;
	//	return 0;
      }
    }
  }
  Gx = Gx_temp; Gz = Gz_temp; Cx = Cx_temp; Cz = Cz_temp;

  //  if ( flag_find_good_code){
  if ( debug ) std::cout<<"Gx 1st row:"<<Gx.get_row(0)<<std::endl; // for debug the random seed
  //make sure both Gx and Gz are full rank
  if ( Gx.row_rank() < Gx.rows() ) {
    if (debug ) std::cout<<"getGoodQuantumCode: Gx not full rank. now make it full rank"<<std::endl;
    Gx = common::make_it_full_rank(Gx);
  }  
  if ( Gz.row_rank() < Gz.rows() ) {
    if (debug) std::cout<<"getGoodQuantumCode: Gz not full rank. now make it full rank"<<std::endl;
    Gz = common::make_it_full_rank(Gz);
  }

  if ( debug) if ( ! flag_find_good_code ) std::cout<<common::color_text("didn't find good code after ")<<search_trial<<" trials"<<std::endl;


  return 0;
}

void set_submatrix(itpp::GF2mat & G, itpp::GF2mat sub, int row, int col){
  //put sub into G, start from (row,col)
  for ( int i =0 ; i < sub.rows(); i ++)
    for ( int j = 0; j< sub.cols(); j++ ){
      G.set(i+row, j+col, sub.get(i,j));
    }
  return;
}


int is_row_reduced_echelon_form(itpp::GF2mat & alpha_Gaz, int debug = 0){
  //check it column by column, from bottom to top
  int get_one=0; //flag on if hit one in that column
  int position_one=-1;//position for one in that column
  int columns_one=0;//columns has a single one. exit for loop when reach alpha_Gaz.rows()
  for ( int i = 0; i<alpha_Gaz.cols();i++){
    get_one=0;
    for ( int j = alpha_Gaz.rows()-1; j > -1; j--){
      if (alpha_Gaz.get(j,i)){
	if (get_one){
	  //get one twice in that column
	  if (debug) std::cout<<"get one twice in that column i="<<i<<std::endl;
	  //	  std::cout<<"*";
	  return 0;
	}else{
	  if (j <= position_one){
	    //skip this column, this column is not independent
	    if (debug) std::cout<<"break the inner for loop for dependent column i = "<<i<<std::endl;
	    break;
	    //	    continue;
	  }else {
	    get_one=1;
	    position_one=j;
	    columns_one++;
	  }
	}
      }
    }
    //if (debug) std::cout<<"broke the inner for loop"<<std::endl;
    if ( columns_one == alpha_Gaz.rows() ){
      break;
    }
  }
  if ( columns_one < alpha_Gaz.rows() ){
    if (debug) std::cout<<"columns_one:"<<columns_one<<" is not full rank"<<alpha_Gaz.rows()<<std::endl;
    //    std::cout<<"*";
    return 0;
  }
  return 1;
}

// generate all code with size na systematically
int generate_code(itpp::GF2mat & Gax, itpp::GF2mat & Gaz, int na, int Gax_row, int id_Gax, int Gaz_row, int id_Gaz, int debug){
  if (debug) std::cout<<na<<","<<Gax_row<<","<<Gaz_row<<std::endl;
  //sanity check
  if (Gaz_row+Gax_row > na-1){
    std::cout<<"generate_code: no logical qubit"<<std::endl;
    throw 2;
  }
  const int id_Gax_MAX = (int) pow(2,  Gax_row * (na-Gax_row) ) -1 ; //maximun all one
  if ( id_Gax <1 || id_Gax > id_Gax_MAX ) {
    std::cout<<"illegal id_Gax: "<<id_Gax<<", id_Gax_MAX = "<<id_Gax_MAX<<std::endl;
    throw 2;
  }
  const int id_Gaz_MAX = (int) pow(2, Gaz_row*(na - Gax_row)) - 1; //maximun all one
  if ( id_Gaz < 1 || id_Gaz > id_Gaz_MAX ){
    std::cout<<"illegal id_Gaz: "<<id_Gaz<<", id_Gaz_MAX = "<<id_Gaz_MAX<<std::endl;
    throw 2;
  }
  //remove duplicate cases for id_Gax and id_Gaz

  itpp::GF2mat beta_Gaz = itpp::GF2mat(itpp::dec2bin(Gaz_row*(na-Gax_row),id_Gaz),false);
  itpp::GF2mat alpha_Gaz(Gaz_row, na-Gax_row);
  if ( debug ) std::cout<<"beta_Gaz = "<<beta_Gaz<<std::endl;
  for ( int i =0;i<Gaz_row;i++){
    if (debug) std::cout<<"set submatrix i = "<<i<<std::endl<<beta_Gaz.get_submatrix(0,i*(na-Gax_row),0, (i+1)*(na-Gax_row)-1)<<std::endl;
    set_submatrix(alpha_Gaz, beta_Gaz.get_submatrix(0,i*(na-Gax_row),0, (i+1)*(na-Gax_row)-1), i,0);
  } 
  if (debug) std::cout<<"alpha_Gaz"<<alpha_Gaz<<std::endl;
  /* decreasing order is ensured in reduced row echelon form, hence not checked here
  for ( int i =0;i<Gaz_row-1;i++){
    if ( itpp::bin2dec(alpha_Gaz.get_row(i)) <= itpp::bin2dec(alpha_Gaz.get_row(i+1))){
      if (debug) std::cout<< "duplicate Gaz with this id_Gaz. no calculation needed. alpha_Gaz/id_Gaz must be in decreasing order"<<std::endl;
      return 2;
    }
    }*/
  //make sure alpa_Gaz is in reduce row echelon form, to remove duplicate cases. return 2 if not in the form
  //this duplicate the check to make sure alpha_Gaz is in decreasing order
  if ( ! is_row_reduced_echelon_form( alpha_Gaz, debug) ) return 2;


  //finish check
 


  Gax = itpp::GF2mat(Gax_row,na);
  // identity matrix in the left part to make it reduce row echelon form.
  set_submatrix(Gax,itpp::gf2dense_eye(Gax_row),0,0);
  //  if (debug) std::cout<<"Gax"<<Gax<<std::endl;


  itpp::GF2mat alpha_Gax = itpp::GF2mat( itpp::dec2bin(Gax_row*(na-Gax_row), id_Gax), false);//false for row vector
  if (debug) std::cout<<"alpha_Gax give the right part of Gax"<<std::endl<<alpha_Gax<<std::endl;
  for ( int i = 0 ; i < Gax_row; i++){
    set_submatrix(Gax,alpha_Gax.get_submatrix(0, i*(na-Gax_row), 0, (i+1)*(na-Gax_row)-1), i, Gax_row);
  }
  if (debug) std::cout<<"Gax"<<Gax<<std::endl;


  //remove duplicate in id_Gax. They could be equal, but permute any two rows give equivalent code, so enfore all rows ( in the right part ) in decreasing order
  for ( int i =0;i<Gax_row-1;i++){
    if ( itpp::bin2dec(Gax.get_submatrix(0,Gax_row,Gax_row-1,na-1).get_row(i)) 
	 < itpp::bin2dec(Gax.get_submatrix(0,Gax_row,Gax_row-1,na-1).get_row(i+1)) ){
      //itpp::bin2dec(alpha_Gaz.get_row(i+1))){
      if (debug) std::cout<< "duplicate Gax with this id_Gax. no calculation needed. id_Gax must be in decreasing order. zero allowed"<<std::endl;
      return 2;
    }
  }
  //check singleton in Gax: row weight = 1
  itpp::bvec bvec_zero = itpp::zeros_b(na);
  for ( int i = 0; i < Gax_row; i++){
    if ( itpp::BERC::count_errors(bvec_zero, Gax.get_row(i)) == 1){
      //      std::cout<<".";
      return 2;
    }
  }

  itpp::GF2mat H = common::nullSpace(Gax);
  if (debug) std::cout<<"nullSpace: H"<<H<<std::endl;

  //check singleton in H: row weight = 1 
  for ( int i = 0; i < na - Gax_row; i++){
    if ( itpp::BERC::count_errors(bvec_zero, H.get_row(i)) == 1){
      //      std::cout<<"+";
      return 2;
    }
  }




  //check id_Gaz


  //  if (debug) std::cout<<"rows_to_remove: "<<rows_to_remove<<std::endl;
  //  remove_rows(&H, rows_to_remove );
  if (debug) std::cout<<"alpha_Gaz"<<alpha_Gaz<<std::endl;
  Gaz=alpha_Gaz*H;
    //itpp::GF2mat(itpp::dec2bin(), false);


  //check singleton in Gaz: col weight = 0
  //  itpp::GF2mat Haz = nullSpace(Gaz);
  itpp::bvec bvec_zero_col=itpp::zeros_b(Gaz.rows());
  for ( int i = 0; i < Gaz.cols(); i++){
    if ( itpp::BERC::count_errors(bvec_zero_col, Gaz.get_col(i)) == 0){
      //      std::cout<<"+";
      return 2;
    }
  }
  
  //  Gaz = alpha_Gaz*H;
  if (debug) std::cout<<"Gaz"<<Gaz<<std::endl;
  return 0;
}


int generate_code(CSSCode & code, int debug){
  return generate_code(code.Gx, code.Gz, code.n, code.Gx_row, code.id_Gx, code.Gz_row, code.id_Gz, debug);
  }




// a version include both reduce and concatenation
// mode=1 for reduce/subsystem product
// mode=2 for concatenation
//only dz is checked cause dx is known to be tight
int product(itpp::GF2mat Gax, itpp::GF2mat Gaz, itpp::GF2mat Gbx, itpp::GF2mat Gbz,int ddax,int ddaz,int ddbx,int ddbz, int debug, int mode){
  //construct code C and calculate the distance; Compare it with the input (estimated) value
  int na=Gax.cols(),nb=Gbx.cols();//,nc=na*nb;//size of the codes

  itpp::GF2mat Cax = common::getC(Gax,Gaz),  Cbx = common::getC(Gbx,Gbz);//This line doesn't allow C to be empty
  itpp::GF2mat Caz = common::getC(Gax,Gaz,1),Cbz = common::getC(Gbx,Gbz,1);//This line doesn't allow C to be empty  
  
  //  Gcz=make_it_full_rank(Gcz);//not needed for calculating distance
  itpp::GF2mat Gcx,Gcz;

  switch ( mode ){
  case 0://reduce/subsystem product, x distance
    
    Gcx = common::kron(Gax,itpp::gf2dense_eye(nb)).concatenate_vertical(common::kron(itpp::gf2dense_eye(na),Gbx));
    Gcz=common::kron(Gaz,Gbz).concatenate_vertical(
					   common::kron(Caz,Gbz)
					   .concatenate_vertical(common::kron(Gaz,Cbz))
					   );
    break;
  case 1://reduce/subsystem product, z distance
    Gcz = common::kron(Gaz,itpp::gf2dense_eye(nb)).concatenate_vertical(common::kron(itpp::gf2dense_eye(na),Gbz));
    Gcx=common::kron(Gax,Gbx)
      .concatenate_vertical(
			    common::kron(Cax,Gbx)
			    .concatenate_vertical(common::kron(Gax,Cbx))
			    );
    break;
  case 2://concatenation
    Gcz = common::kron(Gaz,Cbz).concatenate_vertical(common::kron(itpp::gf2dense_eye(na),Gbz));
    Gcx=common::kron(itpp::gf2dense_eye(na),Gbx).concatenate_vertical( common::kron(Gax,Cbx)   );
    break;
  case 3:
    // chain complex to two CSS codes.
  case 4:
    {
    // chain complex to two CSS codes.
    Gcx=common::kron(Gaz.transpose(), itpp::gf2dense_eye(Gbx.rows()))
      .concatenate_horizontal(common::kron(itpp::gf2dense_eye(Gax.cols()),Gbx))
      .concatenate_horizontal(common::kron(itpp::GF2mat(Gax.cols(),Gax.rows()),itpp::GF2mat(Gbx.rows(),Gbz.rows())));
    Gcx = Gcx
      .concatenate_vertical(
			    common::kron(itpp::GF2mat(Gax.rows(),Gaz.rows()),itpp::GF2mat(Gbx.cols(),Gbx.rows()))
			    .concatenate_horizontal(common::kron(Gax, itpp::gf2dense_eye(Gbx.cols())))
			    .concatenate_horizontal(common::kron(itpp::gf2dense_eye(Gax.rows()),Gbz.transpose()))
			    );
    Gcz=common::kron(itpp::gf2dense_eye(Gaz.rows()), Gbx.transpose())
      .concatenate_horizontal(common::kron(Gaz,itpp::gf2dense_eye(Gbx.cols())))
      .concatenate_horizontal(common::kron(itpp::GF2mat(Gaz.rows(),Gax.rows()),itpp::GF2mat(Gbz.cols(),Gbz.rows())));
    Gcz=Gcz
      .concatenate_vertical(
			    common::kron(itpp::GF2mat( Gaz.cols(),Gaz.rows() ), itpp::GF2mat( Gbz.rows(),Gbx.rows() ))
			    .concatenate_horizontal(common::kron(itpp::gf2dense_eye(Gaz.cols()),Gbz))
			    .concatenate_horizontal(common::kron(Gax.transpose(),itpp::gf2dense_eye(Gbz.rows())))
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
      std::cout<<"mode (3)"<<std::endl;
      break;
    case 4:
      std::cout<<"mode (4)"<<std::endl;
      break;
    }
  }
  if (debug){common::GF2matPrint(Gcx,"Gcx"); common::GF2matPrint(Gcz,"Gcz");}

  if ( ! (Gcx*Gcz.transpose()).is_zero() ){
    std::cout<<"concatenation_lib: not a quantum code "<<std::endl;
    throw "not a quantum code";
  }else{
    if ( debug) std::cout<<"mode ("<<mode<<") is quantum code"<<std::endl;
  }

  //  Gcx=make_it_full_rank(Gcx);//not sure if I need it here
  /*  int daz = quantum_dist(Gax,Gaz,ddaz,1);
  std::cout<<"daz="<<daz<<",ddaz="<<ddaz<<std::endl;
  int dbz = quantum_dist(Gbx,Gbz,ddbz,1);
  std::cout<<"dbz="<<dbz<<",ddbz="<<ddbz<<std::endl;*/
  /*  if (is_quantum_code(Gcx,Gcz) ){
    std::cout<<"C is a quantum Code."<<std::endl;
    }*/
  //  if ( debug ) std::cout<<"Gcx"<<Gcx<<"Gcz"<<Gcz<<std::endl;

  switch ( flag_dist_flip ){
  case 0:
    //x distance
    {
      //  if ( debug ) std::cout<<"Gcx"<<Gcx<<"Gcz"<<Gcz<<std::endl;
      int dax=ddax,dbx=ddbx;
      int dcx = common::quantum_dist(Gcx,Gcz,dax*dbx,debug,0);//donot use estimated value ddax and ddbx
      if (debug) std::cout<<"dax,daz,dbx,dbz = "<<ddax<<","<<ddaz<<","<<ddbx<<","<<ddbz<<","<<std::endl;    
      if (dcx == dax*dbx){
	if (debug) std::cout<<"dcx = dax*dbx = "<<dcx<<std::endl;
	return 0;
      }else if(dcx == common::INF) {
	if (debug) std::cout<<"dcx = "<<dcx<<", dax = "<<dax<<", dbx = "<<dbx<<std::endl;
	return 1;
      }else{
	if (dcx > dax*dbx) std::cout<<"PSEUDO ";
	std::cout<<common::red_text("CASE:")<<" mode ("<<mode<<") dax*dbx="<<dax*dbx<<", dcx="<<dcx;
	std::cout<<". dax,daz,dbx,dbz = "<<ddax<<","<<ddaz<<","<<ddbx<<","<<ddbz<<";";    
	std::cout<<"na,nb,nc,"<<Gax.cols()<<","<<Gbx.cols()<<","<<Gcx.cols()<<";";
	std::cout<<"ka,kb="<<Cax.rows()<<","<<Cbx.rows()<<";";
	return 2;
      }
    }

    break;
  case 1:
    //z distance
    {
      int daz=ddaz,dbz=ddbz;
      int dcz = common::quantum_dist(Gcx,Gcz,daz*dbz,debug,1);//donot use estimated value ddaz and ddbz
      if (debug) std::cout<<"dax,daz,dbx,dbz = "<<ddax<<","<<ddaz<<","<<ddbx<<","<<ddbz<<","<<std::endl;    
      if (dcz == daz*dbz){
	if (debug) std::cout<<"dcz = daz*dbz = "<<dcz<<std::endl;
	return 0;
      }else if(dcz == common::INF) {
	if (debug) std::cout<<"dcz = "<<dcz<<", daz = "<<daz<<", dbz = "<<dbz<<std::endl;
	return 1;
      }else{
	if (dcz > daz*dbz) std::cout<<"PSEUDO ";
	std::cout<<common::red_text("CASE:")<<" mode ("<<mode<<") daz*dbz="<<daz*dbz<<", dcz="<<dcz;
	std::cout<<". dax,daz,dbx,dbz = "<<ddax<<","<<ddaz<<","<<ddbx<<","<<ddbz<<";";
	std::cout<<"na,nb,nc="<<Gax.cols()<<","<<Gbx.cols()<<","<<Gcx.cols()<<";";    
	std::cout<<"ka,kb="<<Cax.rows()<<","<<Cbx.rows()<<";";
	return 2;
      }
    }
  }
  return 0;
  
}



// ============    implementations   ==============


//bool CSSCode::decode(itpp::GF2mat& Gx_temp,itpp::GF2mat& Gz_temp, itpp::bvec e_t, const int perm_try){
bool CSSCode::decode(itpp::bvec e_t, const int perm_try, const int debug){
    //G and S must be full rank
    itpp::GF2mat G=Gz,S=Gx;  
    itpp::GF2mat H=S;    
    H.set_size(S.rows(),S.cols()+1,true);//H=(H_x,s)//H=(H_x,H_z,s)
    G.set_size(G.rows()+1,G.cols(),true);//add one row, use G here because we will use diff=e+X later


    itpp::GF2mat T,U;
    itpp::ivec P;//used for all the gaussian elimination in this file, T_fact(T,U,P)
    itpp::ivec perm(S.cols()+1);
    perm.set(S.cols(),S.cols());//the last col for syndrome, is fixed
    itpp::bvec s=S*e_t;//syndrome s=s_x;//(s_x,s_z);
    H.set_col(H.cols()-1,s);//add syndrome
    //cout<<"parity check matrix H=(S_x,S_z,s). The last column is the syndrome. "<<H<<endl;
    
    int wmin=e_t.length();//minimum weight
    itpp::bvec e_d = itpp::zeros_b(e_t.length() );//error detected with minimum weight wmin //e_d =(e_z,e_x)
    for (int i3=0;i3<perm_try;i3++){//find the error with minimum weight
      //set up random permutation vector for H
      perm.set_subvector(0,itpp::sort_index(itpp::randu( S.cols() )) );

      H.permute_cols( perm, false);
      H.transpose().T_fact(T,U,P);//can I use same T,U,P matrix/vec here despite different size of them: Yes, any size and any value won't affect the result
      //The rows of this matrix are respectively: all errors with non-zero syndrome; all errors with zero syndrome (gauge errors and logical errors); error with the given syndrome. This maybe useful.
      //T should include all the other inequivelant errors. The resize make it cleaner but may lose some meaningful info
      itpp::GF2mat Q=T.get_submatrix(H.rows(),0,H.cols()-1,H.cols()-1);//does permutation matter?  It doesn't matter cause we only look at the zero rows in the H_z.transpose() after gauss. I don't need the nonzero rows to be triangle or not
      //In fact, this Q=(Q_z,Q_x,1) should be Q_tilde. But it is up to the definition, and doesn't affect the permutation
      //cout<<(H*Q.transpose()).density()<<endl;
      //cout<<(H*Q.transpose())<<endl;
      Q.permute_cols(perm,true);//permute back
      /*check permutation relation of Q,H,G
	cout<<"G*Q^T, "<<
	  (G.get_submatrix(0,0,G.rows()-2,G.cols()-1)*
	    Q.get_submatrix(0,0,Q.rows()-2,Q.cols()-2).transpose()
	       ).density()
	       <<endl;
	       cout<<"rank of G: "<<G.get_submatrix(0,0,G.rows()-2,G.cols()-1).T_fact(T,U,P)<<endl;
	       cout<<"rank of Q: "<<Q.get_submatrix(0,0,Q.rows()-2,Q.cols()-2).transpose().T_fact(T,U,P)<<endl;
      */
      
      //cout<<"Q, "<<Q<<endl;
      //cout<<"Q: The last column means if the error match zero syndrome or match the given syndrome in H. Hence only the last row is the error we want."<<endl;

      //get error and check
      itpp::bvec X_t=Q.get_row(Q.rows()-1);//the error detected, X_t=(X_z,X_x,1)
      //      if (debug) std::cout<<"X_t="<<X_t<<std::endl;
      //      if (debug) std::cout<<"Q="<<Q<<std::endl;

      //get the row with minimum weight. add X_t to make sure the last element is one
      int ww=weight(X_t);
      for (int a = 0;a<Q.rows()-1;a++){
	if (ww> weight(X_t+Q.get_row(a)) ){
	  ww= weight(X_t+Q.get_row(a));
	  X_t=X_t+Q.get_row(a);
	}
      }
      X_t.set_size(X_t.size()-1,true);//remove the last 1.  X_t=(X_z,X_x)
      int w=weight(X_t);
      if(wmin>w){//find another error with smaller weight, update it
	wmin=w;
	e_d=X_t;
      }
      //      if (debug) std::cout<<"wmin="<<wmin<<std::endl;
      if (wmin ==0) break;

      H.permute_cols(perm, true);//permute back for another run
    }
  
    //check the error and syndrome    
    itpp::bvec diff_t=e_d+e_t;//same format, (e_z,e_x)
    //    if (debug) std::cout<<"diff_t="<<diff_t<<std::endl;
      
    // bvec diff=get_tilde(diff_t);
    //      G.set_row(G.rows()-1,diff);//check if the diff belong to gauge group

    G.set_row(G.rows()-1,diff_t);//check if the diff belong to gauge group
    //    if (debug) std::cout<<G<<std::endl;
    //    if (debug)  std::cout<<G.rows()<<", rank of new G = "<<G.row_rank()<<std::endl;
    if(G.rows()==G.row_rank()){//not belond to gauge -> belongs to logical group -> bad error
      //      e_bad++;
      return false;
      //cout<<"BAD* ";
      //cout<<"weight(e_t) = "<<weight(e_t)<<", wmin = "<<wmin<<", weight(diff_t) = "<<weight(diff_t)<<endl;
      //      E_input_bad = append_vector(E_input_bad,e_t);
      //      E_output_bad = append_vector(E_output_bad,e_d);
      
    }else{
      return true;
      //      std::cout<<e_d<<", "<<e_t<<", "<<diff_t;
      //      std::cout<<"__good__"<<std::endl;
      //      E_input_good = append_vector(E_input_good,e_t);
      //      E_output_good = append_vector(E_output_good,e_d);
      
    }
    //cout<<S*diff_t<<endl;//check if get the same syndrome
    std::cout<<"should never reach here"<<std::endl;
    return false;
}

//double simulate(itpp::GF2mat Gx, itpp::GF2mat Gz, double p){
double CSSCode::simulate(double p, const int e_try, const int num_cores, const int debug){
  //  itpp::GF2mat Gx,Gz;
  //Gx and Gz must be full rank
  //decode X type error, syndrom s=Gz*e^T

  itpp::RNG_randomize();//get randome seed 
  //#set up parameters#
  //  const int e_try=30000;//number of random errors generated
  const int perm_try=100;//20;//5;//number of trails of random window / permutation;
  int N=Gx.cols();///2;//number of qubits, size of the lattice

  //cout<<"test commutation: G_t*(S.transpoze()) is "<<G_t*(S.transpose())<<endl;
  //define H=(S_x|S_z|s)
  //find Q, the dual of H, which is e+G+L, HQ^T=0

 

  //read error from file
  //  GF2mat E_input=MM_to_GF2mat(filename_E);//here the input is the output of the nonconvergent cases after BP decoding. The first row of this vector is an extra zero vector.
  //int e_try=E_input.rows()-1;//number of total input errors
  //   e_try=(10<e_try) ? 10 :e_try; limit it to 10 to smaller
  //GF2mat E_input_good(e_t,false),E_input_bad(e_t,false),E_output_good(e_t,false),E_output_bad(e_t,false);//false for row vectors
  

  //  if (debug) std::cout<<"before omp, num_cores="<<num_cores<<std::endl;
  int e_bad=0;//count of bad errors
//add pragma here for e_try
//  const int num_cores=32;

#pragma omp parallel for schedule(guided) num_threads(num_cores)
  for(int i1=0;i1<e_try;i1++){
    itpp::bvec e_t = itpp::zeros_b(N);//e_t(2*N);//e_tilde=e_z//(e_z,e_x)
    /*the bug is here
     * for (int i2=0;i2<2*N;i2++)
     *  should be i2<N
    */
    for (int i2=0;i2<N;i2++){//setup random error with error rate p
      e_t.set(i2,(itpp::randu()-p<0)? 1:0); 
    }

    bool decode_result = false;
    if (weight(e_t) == 0){//no need to decode for zero error
      decode_result = true;      
    }else{
      decode_result = decode(e_t, perm_try, debug);
    }

#pragma omp critical
    {
      if (! decode_result){
	e_bad ++;
      }
    }//#pragma omp critical
    //counting e_bad here.
  }//#pragma omp parallel for

  double failure_rate=1.0*e_bad/e_try;

  //print result
  std::cout<<"N = "<<N
      <<", e_try= "<<e_try
      <<", p = "<<p
      <<", failure_rate = "<<failure_rate
	   <<std::endl;
    
  return failure_rate;
}



