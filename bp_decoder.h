//Weilei Zeng, April 8 , 2020
//copied from bp.c

//#include "weilei_lib/my_lib.h"
#include "dist.h"
#include <itpp/itbase.h>
#include <itpp/itcomm.h>
#include <stdio.h>
#include "weilei_lib.h"
//#include <cmath>
#include <stdexcept> //for invalid argument
//using namespace itpp;
//using namespace std;

int summation(itpp::bvec u){
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
  double alpha=1.0; // used for normalized decoder, alpha=1.25
  int decode_mode = 2;
  std::string decode_mode_str="min sum";
  int nvar, ncheck, nedge;

  bool is_initialized=false;
  bool is_H_valid(itpp::GF2mat H_temp);
  void init(itpp::GF2mat H_temp);

  void set_exit_iteration(int exit_iteration_temp);
  void set_decode_mode(int decode_mode_temp);
  void set_decode_mode_str(std::string decode_mode_str_temp);
  void set_silent_mode(bool silent_mode_temp);
  void set_debug_mode(bool debug_mode_temp);
  void print_info();

  
  int decode(itpp::bvec syndrome, const itpp::vec & LLRin, itpp::vec &LLRout);
  int bp_syndrome_llr(itpp::bvec syndrome,  const itpp::vec & LLRin, itpp::vec   & LLRout);
  int bp_schedule(itpp::bvec syndrome,  const itpp::vec & LLRin, itpp::vec   & LLRout);
  bool match_syndrome(itpp::vec LLR, itpp::bvec syndrome);
  int bp_flexible( itpp::bvec syndrome,  const itpp::vec & LLRin, itpp::vec   & LLRout);
  
  int schedule_mode=0;
  std::string schedule_mode_str="no schedule, parallel";
  itpp::mat schedule;
  void set_schedule_mode(int schedule_mode_temp);
  //  void set_schedule(int schedule_mode_temp);
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
  nedge=0;
  for ( int i =0; i<ncheck; i++){
    for ( int j = 0; j<nvar;j++){
      if (H(i,j)) nedge++;
    }
  }
  return;
}

void BP_Decoder::print_info(){
  if (is_initialized){
    std::cout<<"--- BP_Decoder --- nvar = "<<nvar
	     <<", ncheck ="<<ncheck
	     <<", schedule_mode ("<<schedule_mode<<") "<<schedule_mode_str.c_str()
	     <<", decode mode ("<<decode_mode<<") "<<decode_mode_str.c_str()
	     <<", exit_iteration = "<<exit_iteration
	     <<std::endl;

  }else{
    std::cout<<"decoder is not initialized"<<std::endl;
  }
  return;

}

void BP_Decoder::set_exit_iteration(int exit_iteration_temp){
  exit_iteration = exit_iteration_temp;
  return;
}




void BP_Decoder::set_decode_mode(int decode_mode_temp){
  decode_mode = decode_mode_temp;
  alpha=1.0;//reset first
  switch (decode_mode_temp){
  case 1:
    decode_mode_str = "standard";
    break;
  case 2:
    decode_mode_str = "min sum";
    break;
  case 3:// part of min sum, just change alpha
    decode_mode_str = "normalization";
    alpha = 1.25;
    break;
  default:
    throw std::invalid_argument( "BP_Decoder: illegal decode mode" );
  }
  std::cout<<"BP_Dcoder: set decode mode "<<decode_mode<<" - "<<decode_mode_str.c_str()<<std::endl;
  return;
}

   
void BP_Decoder::set_decode_mode_str(std::string decode_mode_str_temp){
  decode_mode_str = decode_mode_str_temp;
  if (decode_mode_str_temp =="standard")    decode_mode = 1;   
  else if (decode_mode_str_temp == "min sum") decode_mode = 2;  
  else if (decode_mode_str_temp == "normalization") decode_mode = 3;
  else throw std::invalid_argument( "BP_Decoder: illegal decode mode string" );
  
  if (debug)  std::cout<<"BP_Dcoder: set decode mode "<<decode_mode<<" - "<<decode_mode_str.c_str()<<std::endl;
  return;
}

void BP_Decoder::set_schedule_mode(int schedule_mode_temp){
  schedule_mode = schedule_mode_temp;
  int s;
  //schedule mode is used in decode();
  switch ( schedule_mode ){
  case 0:
    //no schedule, parrallel schedule
    //    set_decode_mode_str("standard");
    schedule_mode_str = "no schedule, parallel; using bp_syndrome_llr()";
    break;
  case 1: // default schedule mode 1, edge by edge
    //set_decode_mode_str("standard");
    schedule_mode_str  = "edge by edge, using bp_schedule()";
    break;
  case 2: //flexible
    {// edge by edge
      schedule_mode_str ="edge by edge, using bp_flexible()";
      schedule.set_size(2*nedge,3);
      schedule.zeros();
      s=0;
      for ( int i =0; i<ncheck; i++ ){
	for ( int j = 0; j<nvar; j ++){
	  if (H(i,j)){
	    // 0 for c->v, 1 for v->c
	    schedule.set(s,0,0);
	    schedule.set(s,1,i);
	    schedule.set(s,2,j);
	    schedule.set(s+1,0,1);
	    schedule.set(s+1,1,i);
	    schedule.set(s+1,2,j);
	    s +=2;
	  }
	}
      }
      break;
    }
  case 3: //flexible
    { // edge by edge, switch var and check
      schedule_mode_str ="edge by edge, switch var and check order. using bp_flexible()";
      schedule.set_size(2*nedge,3);
      schedule.zeros();
      s=0;
      for ( int j = 0; j<nvar; j ++){
	for ( int i =0; i<ncheck; i++ ){	
	  if (H(i,j)){
	    schedule.set(s,0,0);
	    schedule.set(s,1,i);
	    schedule.set(s,2,j);	
	    schedule.set(s+1,0,1);
	    schedule.set(s+1,1,i);
	    schedule.set(s+1,2,j);
	    s +=2;
	  }
	}
      }
      break;
    }
  case 4: //flexible
    { // random
      set_schedule_mode(2);
      schedule_mode = 4;
      schedule_mode_str ="random purterbation of of mode 2, using bp_flexible()";
      //random permutate
      int permutation=2*nedge;//number of swaps
      itpp::RNG_randomize();//get randome seed 
      itpp::ivec source=itpp::randi(permutation,0,nedge*2-1);
      itpp::ivec target=itpp::randi(permutation,0,nedge*2-1);
      for ( int i =0; i< permutation; i++){
	schedule.swap_rows(source(i),target(i));
      }
      break;
    }
  case 5: //flexible
    { // edge by edge, 
      schedule_mode_str ="edge by edge, follow the paper. using bp_flexible()";
      schedule.set_size(2*nedge,3);
      schedule.zeros();
      s=0;
      std::cout<<"start updating schedule"<<std::endl;
      for ( int j = 0; j<nvar; j ++){
	for ( int i =0; i<ncheck; i++ ){	
	  if (H(i,j)){
	    schedule.set(s,0,0);
	    schedule.set(s,1,i);
	    schedule.set(s,2,j);
	    s ++;
	  }	
	}
	for ( int i =0; i<ncheck; i++ ){	
	  if (H(i,j)){
	    schedule.set(s,0,1);
	    schedule.set(s,1,i);
	    schedule.set(s,2,j);
	    s ++;
	  }	
	}
      }
      std::cout<<"finish updating schedule"<<std::endl;
      break;
    }
    
  default:
    throw std::invalid_argument( "BP_Decoder: illegal schedule" );
  }
  return;
}


int BP_Decoder::decode( itpp::bvec syndrome,  const itpp::vec & LLRin, itpp::vec   & LLRout){
  // a wrapper to determine using which decoding function
  switch ( schedule_mode){
  case 0: // no schedule
    return bp_syndrome_llr(syndrome,  LLRin,  LLRout);
  case 1://same position for u and v, one by one
    return bp_schedule( syndrome, LLRin, LLRout);
  case 2://same position for u and v, one by one
  case 3://same position for v and u, one by one,
  case 4: //random permutation of case 2
  case 5: //follow the paper
    return bp_flexible( syndrome, LLRin, LLRout);
  default:
    throw std::invalid_argument( "BP_Decoder: illegal schedule" );
    //return bp_syndrome_llr(syndrome,  LLRin,  LLRout);
  }  
}

bool BP_Decoder::match_syndrome(itpp::vec LLR, itpp::bvec syndrome){
  itpp::bvec  error = LLR < 0;
  return itpp::GF2mat(H*error - syndrome).is_zero();
}

//int BP_Decoder::bp_syndrome_llr(const GF2mat H,  itpp::bvec syndrome,  itpp::vec & LLRin, itpp::vec   & LLRout, int exit_iteration, int decode_mode){
int BP_Decoder::bp_syndrome_llr( itpp::bvec syndrome,  const itpp::vec & LLRin, itpp::vec   & LLRout){  
  // input: syndrome vector s, loglikelihood ratio LLRin and LLRout
  // exit_iteration: max number of iteration
  // decode_mode 1: standard, 2: min sum
  //LLR(x) = log( p(x)/ (1-p(x)) )
  // initially we assume zero error, so LLR = LLR(0)=log ( (1-p)/p )>0
  // bits_out = LLRout < 0;
  //output: number of iteration, negative if not converge.

  if ( itpp::GF2mat(syndrome).is_zero() ){
    //return zero error vector, which is the default input for LLRout
    return 1;
  }
  
  // bool debug = false;// enable/disable printing
  //initialize
  //int nvar = H.cols(), ncheck = H.rows();
  if (debug) std::cout<<"nvar = "<<nvar <<", ncheck = "<<ncheck<<std::endl;
  itpp::mat llrs = itpp::zeros(ncheck, nvar), LLRs = itpp::zeros(ncheck, nvar);  
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
  //  if (debug) std::cout<<"finish initialize"<<std::endl;

  int update_count=0;
  double sum=0;
  double llr,LLR;
  std::string str="";
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
		  str += std::to_string(i)+",";
		  str += std::to_string(j)+",";
		  str += std::to_string(k)+",";
		  str += std::to_string(prod)+",";
		  str += std::to_string(llrs(i,k))+",";
		  //if (debug) std::cout<<",prod = "<<prod<<" i="<<i <<std::endl;
		}
	      }
	    }
	    
	    LLR = atanh(prod)*2;
	    if ( syndrome(i) ){
	      LLR = -LLR;
	    }
	    LLRs.set(i,j,LLR);

	    if (debug) if ( std::abs(LLR) > 1000000.0) std::cout<<"LLRs: LLR = "<<LLR<<", prod = "<<prod<<"\n"<<str<<std::endl<<"H.get_row(i)="<<H.get_row(i)<<std::endl <<"llrs.get_row(i)="<<llrs.get_row(i)<<std::endl;
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
		  prod = (llr < prod) ? llr: prod;//min( prod, llr);
			      
		    //		  prod = prod * tanh( llrs(i,k)/2 );
		    /*str += std::to_string(i)+",";
		  str += std::to_string(j)+",";
		  str += std::to_string(k)+",";
		  str += std::to_string(prod)+",";
		  str += std::to_string(llrs(i,k))+",";*/
		  //if (debug) std::cout<<",prod = "<<prod<<" i="<<i <<std::endl;
		}
	      }
	    }
	    LLR = sign * prod;
	    //	    double LLR = atanh(prod)*2;
	    if ( syndrome(i) ){
	      LLR = -LLR;
	    }
	    LLRs.set(i,j,LLR);

	    if (debug) if ( std::abs(LLR) > INF_BP ) std::cout<<"LLRs: LLR = "<<LLR<<", prod = "<<prod<<"\n"<<str<<std::endl<<"H.get_row(i)="<<H.get_row(i)<<std::endl <<"llrs.get_row(i)="<<llrs.get_row(i)<<std::endl;
	  }
	}
      }
    }
    
    //    if (debug) std::cout<<"finish check to variable update"<<std::endl;    
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
	  //	  if ( std::abs(sum) > 1000)	  std::cout<<"llrs: sum = "<<sum<<"\n"<<LLRs.get_col(j)<<std::endl;
	}

      }
    }

    //    if (debug) std::cout<<"finish variable to checkupdate"<<std::endl;
  
    // get output LLRout and check result
    //    match_syndrome = true;
    for ( int j=0; j<nvar; j++){
        sum=LLRin(j);
	//if (debug) std::cout<<" sum = "<<sum<<std::endl;
	for ( int t=0; t<ncheck; t++){
	  //if (debug) std::cout<<"t = "<<t<<", sum = "<<sum<<std::endl;
	  if ( H(t,j) ){
	      sum += LLRs(t,j);
	  }
	}
	//if (debug) std::cout<<"LLRout = "<<LLRout<<std::endl;
	LLRout.set(j,sum);
    }
    if (debug) std::cout<<"update_count = "<<update_count<<", LLRout = "<<floor(LLRout)<<std::endl ;
    //if (debug) std::cout<<"update_count = "<<update_count<<std::endl;
    //if (debug) draw_toric_x_error(LLRout<0);
    update_count++;
    if ( match_syndrome( LLRout, syndrome) ){
      break;
    }        
  }
    
  if (debug) std::cout<<"LLRout = "<<LLRout<<std::endl;

  //  if (debug) std::cout<<"llrs = "<<llrs<<std::endl;
  //if (debug) std::cout<<"LLRs = "<<LLRs<<std::endl;

  //not converge, output negative value
  if (! match_syndrome( LLRout, syndrome) ){
    update_count = - update_count;
  }
    
  return update_count ;
}


int BP_Decoder::bp_schedule( itpp::bvec syndrome,  const itpp::vec & LLRin, itpp::vec   & LLRout){
 // input: syndrome vector s, loglikelihood ratio LLRin and LLRout
  // exit_iteration: max number of iteration
  // decode_mode 1: standard, 2: min sum
  //LLR(x) = log( p(x)/ (1-p(x)) )
  // initially we assume zero error, so LLR = LLR(0)=log ( (1-p)/p )>0
  // bits_out = LLRout < 0;
  //output: number of iteration, negative if not converge.

  if ( itpp::GF2mat(syndrome).is_zero() ){
    //return zero error vector, which is the default input for LLRout
    return 1;
  }
  
  // bool debug = false;// enable/disable printing
  //initialize
  //int nvar = H.cols(), ncheck = H.rows();
  if (debug) std::cout<<"nvar = "<<nvar <<", ncheck = "<<ncheck<<std::endl;
  itpp::mat llrs = itpp::zeros(ncheck, nvar), LLRs = itpp::zeros(ncheck, nvar);  
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
  //  if (debug) std::cout<<"finish initialize"<<std::endl;

  // *********************************
  int update_count=0;
  double sum=0;
  double LLR;//llr is not used yet
  std::string str="";
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
		  /*str += std::to_string(i)+",";
		  str += std::to_string(j)+",";
		  str += std::to_string(k)+",";
		  str += std::to_string(prod)+",";
		  str += std::to_string(llrs(i,k))+",";*/
		  //if (debug) std::cout<<",prod = "<<prod<<" i="<<i <<std::endl;
		}
	      }
	    }
	    
	    LLR = atanh(prod)*2;
	    if ( syndrome(i) ){
	      LLR = -LLR;
	    }
	    LLRs.set(i,j,LLR);

	    if (debug) if ( std::abs(LLR) > 1000000.0) std::cout<<"LLRs: LLR = "<<LLR<<", prod = "<<prod<<"\n"<<str<<std::endl<<"H.get_row(i)="<<H.get_row(i)<<std::endl <<"llrs.get_row(i)="<<llrs.get_row(i)<<std::endl;
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
	  //	  if ( std::abs(sum) > 1000)	  std::cout<<"llrs: sum = "<<sum<<"\n"<<LLRs.get_col(j)<<std::endl;
	}

      }
    }

    //    if (debug) std::cout<<"finish variable to checkupdate"<<std::endl;
  
    // get output LLRout and check result
    //    match_syndrome = true;
    for ( int j=0; j<nvar; j++){
        sum=LLRin(j);
	//if (debug) std::cout<<" sum = "<<sum<<std::endl;
	for ( int t=0; t<ncheck; t++){
	  //if (debug) std::cout<<"t = "<<t<<", sum = "<<sum<<std::endl;
	  if ( H(t,j) ){
	      sum += LLRs(t,j);
	  }
	}
	//if (debug) std::cout<<"LLRout = "<<LLRout<<std::endl;
	LLRout.set(j,sum);
    }
    if (debug) std::cout<<"update_count = "<<update_count<<", LLRout = "<<floor(LLRout)<<std::endl ;
    //if (debug) std::cout<<"update_count = "<<update_count<<std::endl;
    //if (debug) draw_toric_x_error(LLRout<0);
    update_count++;
    if ( match_syndrome( LLRout, syndrome) ){
      break;
    }        
  }
    
  if (debug) std::cout<<"LLRout = "<<LLRout<<std::endl;

  //  if (debug) std::cout<<"llrs = "<<llrs<<std::endl;
  //if (debug) std::cout<<"LLRs = "<<LLRs<<std::endl;

  //not converge, output negative value
  if (! match_syndrome( LLRout, syndrome) ){
    update_count = - update_count;
  }
    
  return update_count ;

}


int BP_Decoder::bp_flexible( itpp::bvec syndrome,  const itpp::vec & LLRin, itpp::vec   & LLRout){
  // a flexible function to encorporate all variations
 // input: syndrome vector s, loglikelihood ratio LLRin and LLRout
  // exit_iteration: max number of iteration
  // decode_mode 1: standard, 2: min sum
  //LLR(x) = log( p(x)/ (1-p(x)) )
  // initially we assume zero error, so LLR = LLR(0)=log ( (1-p)/p )>0
  // bits_out = LLRout < 0;
  //output: number of iteration, negative if not converge.

  if ( itpp::GF2mat(syndrome).is_zero() )     return 1;
    //return zero error vector, which is the default input for LLRout
  
  //initialize
  if (debug) std::cout<<"nvar = "<<nvar <<", ncheck = "<<ncheck<<std::endl;
  itpp::mat llrs = itpp::zeros(ncheck, nvar), LLRs = itpp::zeros(ncheck, nvar);  
  //LLRout.set_size(nvar);should be the same size
  
  for ( int i = 0; i< ncheck ; i++){
    for ( int j=0; j<nvar; j++){
      if (H(i,j)) {
	llrs.set(i,j,LLRin(j));
	//llrs.set(i,j,log( (1-p)/p ));
      }
    }
  }

  // ********************************* start updating cycle
  int update_count=0;
  double sum=0;
  double LLR, llr;
  std::string str="";
  int sign; 
  double prod=1.0;
  while ( update_count < exit_iteration ){
    int i,j;
    for ( int is = 0; is < nedge*2; is ++){
      i = schedule(is,1);
      j = schedule(is,2);      
      //      for ( int i = 0; i< ncheck ; i++){
      //for ( int j=0; j<nvar; j++){
      int direction = schedule(is,0);
      switch ( direction ){
      case 0:     //check to variable update, LLR
	{
	  switch ( decode_mode ) {
	  case 1: //standard
	    {
	      prod=1.0;
	      for ( int k=0; k<nvar; k++){
		if ( H(i,k) ){
		  if ( k != j ) prod = prod * tanh( llrs(i,k)/2 );		  
		}
	      }
	      
	      LLR = atanh(prod)*2;
	      if ( syndrome(i) ) LLR = -LLR;	      
	      LLRs.set(i,j,LLR);
	    
	      if (debug) if ( std::abs(LLR) > 1000000.0) std::cout<<"LLRs: LLR = "<<LLR<<", prod = "<<prod<<"\n"<<str<<std::endl<<"H.get_row(i)="<<H.get_row(i)<<std::endl <<"llrs.get_row(i)="<<llrs.get_row(i)<<std::endl;
	      break;
	    }

	  case 2://min sum
	    //      for ( int i = 0; i< ncheck ; i++){
	    //for ( int j=0; j<nvar; j++){
	    //  if (H(i,j)) {
	    {
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
		    prod = (llr<prod) ? llr:prod;//min( prod, llr);
		  }
		}
	      }
	      LLR = sign * prod;
	      if ( syndrome(i) ) LLR = -LLR;
	      LLRs.set(i,j,LLR);

	      if (debug) if ( std::abs(LLR) > INF_BP ) std::cout<<"LLRs: LLR = "<<LLR<<", prod = "<<prod<<"\n"<<str<<std::endl<<"H.get_row(i)="<<H.get_row(i)<<std::endl <<"llrs.get_row(i)="<<llrs.get_row(i)<<std::endl;

	      break;
	    }
	    break;
	  }
	}
      case 1:
	{
	    // variable to check, llr updating
	  sum=  LLRin(j);		
	  for ( int t=0; t<ncheck; t++){
	    if ( H(t,j) ){
	      if ( t != i ) {
		sum += LLRs(t,j);		
	      }
	    }
	  }
	  llrs.set(i,j,sum);
	  break;
	}
      }
    }
  
    // get output LLRout and check result
    for ( int j=0; j<nvar; j++){
        sum=LLRin(j);
	//if (debug) std::cout<<" sum = "<<sum<<std::endl;
	for ( int t=0; t<ncheck; t++){
	  //if (debug) std::cout<<"t = "<<t<<", sum = "<<sum<<std::endl;
	  if ( H(t,j) ){
	      sum += LLRs(t,j);
	  }
	}
	//if (debug) cout<<"LLRout = "<<LLRout<<std::endl;
	LLRout.set(j,sum);
    }
    if (debug) std::cout<<"update_count = "<<update_count<<", LLRout = "<<floor(LLRout)<<std::endl ;
    //if (debug) std::cout<<"update_count = "<<update_count<<std::endl;
    //if (debug) draw_toric_x_error(LLRout<0);
    update_count++;
    if ( match_syndrome( LLRout, syndrome) ){
      break;
    }        
  }
    
  if (debug) std::cout<<"LLRout = "<<LLRout<<std::endl;

  //  if (debug) std::cout<<"llrs = "<<llrs<<std::endl;
  //if (debug) std::cout<<"LLRs = "<<LLRs<<std::std::endl;

  //not converge, output negative value
  if (! match_syndrome( LLRout, syndrome) ){
    update_count = - update_count;
  }
    
  return update_count ;

}
