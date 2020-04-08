//Weilei Zeng, April 8 , 2020
//copied from bp.c

//#include "weilei_lib/my_lib.h"
#include "dist.h"
#include <itpp/itbase.h>
#include <itpp/itcomm.h>
#include <stdio.h>
#include "my_lib.h"
#include <cmath>
//using namespace itpp;
//using namespace std;

class BP_Decoder{
 public:
  const double INF_BP=10000;
  itpp::GF2mat H;//parity check matrix
  int exit_iteration=50;
  int decode_mode = 1;
  std::string decode_mode_str="standard";
  int nvar, ncheck;
  void init(itpp::GF2mat H_parity_check);
  
  void decode();

  int bp_syndrome_llr(bvec syndrome,  vec & LLRin, vec   & LLRout);
  bool match_syndrome(vec LLR, bvec syndrome);
  
};

void BP_Decoder::init(itpp::GF2mat H_parity_check){
  H = H_parity_check;
  nvar = H.cols();
  ncheck = H.rows();
  return;
}



bool BP_Decoder::match_syndrome(vec LLR, bvec syndrome){
  bvec  error = LLR < 0;
  return GF2mat(H*error - syndrome).is_zero();
}


//int BP_Decoder::bp_syndrome_llr(const GF2mat H,  bvec syndrome,  vec & LLRin, vec   & LLRout, int exit_iteration, int decode_mode){
int BP_Decoder::bp_syndrome_llr( bvec syndrome,  vec & LLRin, vec   & LLRout){  
  // input: parity check matrix H, syndrome vector s, loglikelihood ratio LLRin and LLRout
  // exit_iteration: max number of iteration
  // decode_mode 1: standard, 2: min sum
  //LLR(x) = log( p(x)/ (1-p(x)) )
  // initially we assume zero error, so LLR = LLR(0)=log ( (1-p)/p )>0
  // bits_out = LLRout < 0;
  //output: number of iteration, negative if not converge.

  if ( GF2mat(syndrome).is_zero() ){
    //return zero error vector, which is the default input for LLRout
    return 1;
  }
  
  bool debug = false;// enable/disable printing
  //initialize
  int nvar = H.cols(), ncheck = H.rows();
  if (debug) cout<<"nvar = "<<nvar <<", ncheck = "<<ncheck<<endl;
  mat llrs = zeros(ncheck, nvar), LLRs=zeros(ncheck, nvar);  
  LLRout.set_size(nvar);
  //  bool match_syndrome = false;// a flag to indicate whether the syndrome has been satisfied
  
  for ( int i = 0; i< ncheck ; i++){
    for ( int j=0; j<nvar; j++){
      if (H(i,j)) {
	llrs.set(i,j,LLRin(j));
	//llrs.set(i,j,log( (1-p)/p ));
      }
    }
  }
  //  if (debug) cout<<"finish initialize"<<endl;

  //*********************************
  int update_count=0;
  double sum=0;
  double llr,LLR;
  string str="";
  int sign = 1;
  double prod=1.0;
  //  bool degree_one=true;
  while ( update_count < exit_iteration ){
    //check to variable update, LLR
    switch (decode_mode){
    case 1://standard
      for ( int i = 0; i< ncheck ; i++){
	for ( int j=0; j<nvar; j++){
	  if (H(i,j)) {
	    prod=1.0;
	    str = "prod list: i,j,k,prod,llr:";	  
	    for ( int k=0; k<nvar; k++){
	      if ( H(i,k) ){
		if ( k != j ) {
		  prod = prod * tanh( llrs(i,k)/2 );
		  str += to_string(i)+",";
		  str += to_string(j)+",";
		  str += to_string(k)+",";
		  str += to_string(prod)+",";
		  str += to_string(llrs(i,k))+",";
		  //if (debug) cout<<",prod = "<<prod<<" i="<<i <<endl;
		}
	      }
	    }
	    
	    LLR = atanh(prod)*2;
	    if ( syndrome(i) ){
	      LLR = -LLR;
	    }
	    LLRs.set(i,j,LLR);

	    if (debug) if ( std::abs(LLR) > 1000000.0) cout<<"LLRs: LLR = "<<LLR<<", prod = "<<prod<<"\n"<<str<<endl<<"H.get_row(i)="<<H.get_row(i)<<endl <<"llrs.get_row(i)="<<llrs.get_row(i)<<endl;
	  }
	}
      }
    case 2://min sum
      for ( int i = 0; i< ncheck ; i++){
	for ( int j=0; j<nvar; j++){
	  if (H(i,j)) {
	    prod=INF_BP;
	    sign = 1;
	    //	    str = "prod list: i,j,k,prod,llr:";	  
	    for ( int k=0; k<nvar; k++){
	      if ( H(i,k) ){
		if ( k != j ) {
		  if (llrs(i,k)>0){
		    llr = llrs(i,k);
		  }else{
		    llr = -llrs(i,k);
		    sign = -sign;
		  }
		  prod = min( prod, llr);
			      
		    //		  prod = prod * tanh( llrs(i,k)/2 );
		    /*str += to_string(i)+",";
		  str += to_string(j)+",";
		  str += to_string(k)+",";
		  str += to_string(prod)+",";
		  str += to_string(llrs(i,k))+",";*/
		  //if (debug) cout<<",prod = "<<prod<<" i="<<i <<endl;
		}
	      }
	    }
	    LLR = sign * prod;
	    //	    double LLR = atanh(prod)*2;
	    if ( syndrome(i) ){
	      LLR = -LLR;
	    }
	    LLRs.set(i,j,LLR);

	    if (debug) if ( std::abs(LLR) > 1000000.0) cout<<"LLRs: LLR = "<<LLR<<", prod = "<<prod<<"\n"<<str<<endl<<"H.get_row(i)="<<H.get_row(i)<<endl <<"llrs.get_row(i)="<<llrs.get_row(i)<<endl;
	  }
	}
      }
    }
    
    //    if (debug) cout<<"finish check to variable update"<<endl;
    
    //variable to check update, llr
    


    for ( int i = 0; i< ncheck ; i++){
      for ( int j=0; j<nvar; j++){
	if (H(i,j)) {
	  sum=  LLRin(j);	
	
	  for ( int t=0; t<ncheck; t++){
	    if ( H(t,j) ){
	      if ( t != i ) {
		sum += LLRs(t,j);
		
	      }
	    }
	  }
	  llrs.set(i,j,sum);
	  //	  if ( std::abs(sum) > 1000)	  cout<<"llrs: sum = "<<sum<<"\n"<<LLRs.get_col(j)<<endl;
	}

      }
    }

    //    if (debug) cout<<"finish variable to checkupdate"<<endl;
  
    // get output LLRout and check result
    //    match_syndrome = true;
    for ( int j=0; j<nvar; j++){
        sum=LLRin(j);
	//if (debug) cout<<" sum = "<<sum<<endl;
	for ( int t=0; t<ncheck; t++){
	  //if (debug) cout<<"t = "<<t<<", sum = "<<sum<<endl;
	  if ( H(t,j) ){
	      sum += LLRs(t,j);
	  }
	}
	//if (debug) cout<<"LLRout = "<<LLRout<<endl;
	LLRout.set(j,sum);

	//	if ( std::abs(sum) > 1000)	  cout<<"LLRout: sum = "<<sum<<endl;
	/*if ( std::abs(sum) > 1000){
	  cout<< sum<<endl;
	  }*/	
	//if (debug) cout<<" sum = "<<sum<<endl;
	/*
	if ( sum * (double) syndrome(j) <0 ){ //syndrome not satisfied
	  match_syndrome = false;
	  //break;
	  }*/	
    }
    if (debug) cout<<"update_count = "<<update_count<<", LLRout = "<<floor(LLRout)<<endl ;
    //if (debug) cout<<"update_count = "<<update_count<<endl;
    //if (debug) draw_toric_x_error(LLRout<0);
    update_count++;
    if ( match_syndrome( LLRout, syndrome) ){
      break;
    }        
  }
    
  if (debug) cout<<"LLRout = "<<LLRout<<endl;

  //  if (debug) cout<<"llrs = "<<llrs<<endl;
  //if (debug) cout<<"LLRs = "<<LLRs<<endl;

  //not converge, output negative value
  if (! match_syndrome( LLRout, syndrome) ){
    update_count = - update_count;
  }
    
  return update_count ;
}
