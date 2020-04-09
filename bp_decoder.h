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

int summation(bvec u){
  int result = 0;
  for ( int i=0; i<u.size(); i++){
    if (u(i)) result ++;
  }
  return result;
}


class BP_Decoder{
 public:
  bool silent_mode = false;//not is use
  bool debug = false;
  const double INF_BP=10000;
  itpp::GF2mat H;//parity check matrix
  int exit_iteration=50;
  int decode_mode = 1;
  std::string decode_mode_str="standard";
  int nvar, ncheck;
  bool is_initialized=false;
  bool is_H_valid(itpp::GF2mat H_temp);
  void init(itpp::GF2mat H_temp);
  void set_exit_iteration(int exit_iteration_temp);
  void set_decode_mode(int decode_mode_temp);
  void set_decode_mode_str(std::string decode_mode_str_temp);
  void set_silent_mode(bool silent_mode_temp);
  void set_debug_mode(bool debug_mode_temp);
  void print_info();

  
  int decode(bvec syndrome, const vec & LLRin, vec &LLRout);
  int bp_syndrome_llr(bvec syndrome,  const vec & LLRin, vec   & LLRout);
  int bp_schedule(bvec syndrome,  const vec & LLRin, vec   & LLRout);
  bool match_syndrome(vec LLR, bvec syndrome);

  int schedule_mode=0;
  void set_schedule_mode(int schedule_mode_temp);
  
};

void BP_Decoder::set_silent_mode(bool silent_mode_temp){
  silent_mode = silent_mode_temp;  
  return;
}
void BP_Decoder::set_debug_mode(bool debug_mode_temp){
  debug = debug_mode_temp;
  return;
}


bool BP_Decoder::is_H_valid(itpp::GF2mat H_temp){
  //check degree >=2
  for ( int i =0; i<H_temp.rows();i++){
    if ( summation(H_temp.get_row(i)) < 2 ) return false;      
  }
  //variable degree >=1
  for ( int j = 0; j< H_temp.cols(); j++){
    if ( summation(H_temp.get_col(j)) < 1 ) return false;
  }
  return true;
}

void BP_Decoder::init(itpp::GF2mat H_temp){
  //set up parity check matrix
  H = H_temp;
  nvar = H.cols();
  ncheck = H.rows();  
  if ( ! is_H_valid(H) ) throw std::invalid_argument( "BP_Decoder: invalid parity check matrix H" );
  is_initialized = true;
  return;
}

void BP_Decoder::print_info(){
  if (is_initialized){
    std::cout<<"--- BP_Decoder --- nvar = "<<nvar
	     <<", ncheck ="<<ncheck
	     <<", schedule_mode = "<<schedule_mode
	     <<", decode mode ("<<decode_mode<<") "<<decode_mode_str.c_str()
	     <<", exit_iteration = "<<exit_iteration
	     <<std::endl;

  }else{
    std::cout<<"decoder is not initialized"<<endl;
  }
  return;

}

void BP_Decoder::set_exit_iteration(int exit_iteration_temp){
  exit_iteration = exit_iteration_temp;
  return;
}




void BP_Decoder::set_decode_mode(int decode_mode_temp){
  switch (decode_mode_temp){
  case 1:
    decode_mode = 1;
    decode_mode_str = "standard";
    break;
  case 2:
    decode_mode = 2;
    decode_mode_str = "min sum";
    break;
  default:
    throw std::invalid_argument( "BP_Decoder: illegal decode mode" );
  }
  cout<<"BP_Dcoder: set decode mode "<<decode_mode<<" - "<<decode_mode_str.c_str()<<std::endl;
  return;
}

   
void BP_Decoder::set_decode_mode_str(std::string decode_mode_str_temp){
  if (decode_mode_str_temp =="standard"){
    decode_mode = 1;
    decode_mode_str = "standard";
  }
  else if (decode_mode_str_temp == "min sum"){
    decode_mode = 2;
    decode_mode_str = "min sum";
  }else{  
    throw std::invalid_argument( "BP_Decoder: illegal decode mode string" );
  }
  if (debug)  cout<<"BP_Dcoder: set decode mode "<<decode_mode<<" - "<<decode_mode_str.c_str()<<std::endl;
  return;
}

void BP_Decoder::set_schedule_mode(int schedule_mode_temp){
  schedule_mode = schedule_mode_temp;
  if ( schedule_mode ==1 )  set_decode_mode_str("standard");
}

int BP_Decoder::decode( bvec syndrome,  const vec & LLRin, vec   & LLRout){
  // a wrapper to determine using which decoding function
  switch ( schedule_mode){
  case 1://same position for u and v, one by one
    return bp_schedule( syndrome, LLRin, LLRout);
  default: //no schedule
    return bp_syndrome_llr(syndrome,  LLRin,  LLRout);
  }  
}

bool BP_Decoder::match_syndrome(vec LLR, bvec syndrome){
  bvec  error = LLR < 0;
  return GF2mat(H*error - syndrome).is_zero();
}

//int BP_Decoder::bp_syndrome_llr(const GF2mat H,  bvec syndrome,  vec & LLRin, vec   & LLRout, int exit_iteration, int decode_mode){
int BP_Decoder::bp_syndrome_llr( bvec syndrome,  const vec & LLRin, vec   & LLRout){  
  // input: syndrome vector s, loglikelihood ratio LLRin and LLRout
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
  
  // bool debug = false;// enable/disable printing
  //initialize
  //int nvar = H.cols(), ncheck = H.rows();
  if (debug) cout<<"nvar = "<<nvar <<", ncheck = "<<ncheck<<endl;
  mat llrs = zeros(ncheck, nvar), LLRs=zeros(ncheck, nvar);  
  //LLRout.set_size(nvar);should be the same size
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

	    if (debug) if ( std::abs(LLR) > INF_BP ) cout<<"LLRs: LLR = "<<LLR<<", prod = "<<prod<<"\n"<<str<<endl<<"H.get_row(i)="<<H.get_row(i)<<endl <<"llrs.get_row(i)="<<llrs.get_row(i)<<endl;
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


int BP_Decoder::bp_schedule( bvec syndrome,  const vec & LLRin, vec   & LLRout){
 // input: syndrome vector s, loglikelihood ratio LLRin and LLRout
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
  
  // bool debug = false;// enable/disable printing
  //initialize
  //int nvar = H.cols(), ncheck = H.rows();
  if (debug) cout<<"nvar = "<<nvar <<", ncheck = "<<ncheck<<endl;
  mat llrs = zeros(ncheck, nvar), LLRs=zeros(ncheck, nvar);  
  //LLRout.set_size(nvar);should be the same size
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
  double LLR;//llr is not used yet
  string str="";
  //  int sign; //not used yet
  double prod=1.0;
  //  bool degree_one=true;
  while ( update_count < exit_iteration ){
    //check to variable update, LLR
    //use standard updating rule, no min sum
      for ( int i = 0; i< ncheck ; i++){
	for ( int j=0; j<nvar; j++){
	  if (H(i,j)) {
	    prod=1.0;
	    //	    str = "prod list: i,j,k,prod,llr:";	  
	    for ( int k=0; k<nvar; k++){
	      if ( H(i,k) ){
		if ( k != j ) {
		  prod = prod * tanh( llrs(i,k)/2 );
		  /*str += to_string(i)+",";
		  str += to_string(j)+",";
		  str += to_string(k)+",";
		  str += to_string(prod)+",";
		  str += to_string(llrs(i,k))+",";*/
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
	    //}
	  //}
	//}

    
      //    for ( int i = 0; i< ncheck ; i++){
      //for ( int j=0; j<nvar; j++){
      //if (H(i,j)) {
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
