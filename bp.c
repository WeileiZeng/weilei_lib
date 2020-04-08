//Weilei Zeng, April 28
//return min weight of rows in C, which is the distance of the code
//random window method is applied with default number of permutation 10.

//#include "weilei_lib/my_lib.h"
#include "dist.h"
#include <itpp/itbase.h>
#include <itpp/itcomm.h>
#include <stdio.h>
#include "my_lib.h"
#include <cmath>
using namespace itpp;
using namespace std;



bvec reduce_weight(bvec e, GF2mat G){
  //reduce the weight of en error vector e, G include all cycles (equivalent error). G is assument to be the generator matrix of an Quantum LDPC CSS code. such that rows of G correspond to the smallest cycles. It is assument that e is an effectively small/zero weight error, which contains many trivial cycles. This function will reduce most of those small individual cycles so that we can count the effective weight of e.
  //  complexity of this program: ~ G.rows() * weight(e)
  int row=G.rows();
  bool converge = false;
  bvec e1;
  int wt=weight(e);
  while (! converge){
    converge=true;//assume it is already at min weight
    for (int i =0; i< row;i++){
      e1=e+G.get_row(i);
      if (weight(e1) < wt){
	e=e1;
	wt=weight(e1);
	converge = false;
      }
    }
  }
  return e;
}


bvec qllr_to_bvec(QLLRvec llr, int bound){
  //find the appropriate bound value and return binary vector from this qllr.
  //Its sign is the same as the sign of "bound"
  /*
  int bound1=bound;
  int m  = max( abs(llr));
  int upper = m/100;
  int lower = m/10000;
  bound1=min(upper,bound1);
  bound1=max(lower,bound1);
  bound1 = bound / abs(bound) *bound1;
  
  bound1=0;
  bvec bits_out = llr < bound1;

  cout<<"actual bound1 = "<<bound1<<endl;
  return bits_out;
  */
  bvec bits_out = llr <0;
  return bits_out;
}


//for partial sum project

GF2mat read_one_matrix(char * filename_prefix, char * type, char * filename_suffix){
  
  char filename[255];
  sprintf( filename,"%s%s%s",filename_prefix, type, filename_suffix);//append p in the file name
  //cout<<"read one matrix: "<<filename<<endl;
  return MM_to_GF2mat(filename);
}


int check_matrices(GF2mat * G, GF2mat * H, GF2mat *U, GF2mat * W, mat * K){
  //check if those matrices are valid
  //it seems that H has some weight 1 column, which is not permited in H. bugs find when weight =7, size =9 or 13
  cout<<"debug read matrices by checking some relation of the output matrices"<<endl;
  //  GF2matPrint(*G,"G");
  GF2matPrint(*G,(char *) "G");
  GF2matPrint(*H,(char *) "H");
  GF2matPrint(*U,(char *) "U");
  GF2matPrint(*W,(char *) "W");
  matPrint(*K,(char *) "K");
  
  cout<<"if G*H is zero? ->"<<( ( *G *  (*H).transpose()).is_zero() ?"yes":"no" )<<endl;

  if ( G->cols() == H->cols() && G->cols()==U->cols()  && G->cols() == K->cols()-1 && G->rows() == W->rows() ){
    cout<<"dimension matches"<<endl;
  }else{
    cout<<"dimension does NOT match"<<endl;
  }
  
  
  //check column weight of H, variable to check nodes require min wt >=1
  for ( int i =0; i< H->cols();i++){
    int wt=0;
    for ( int j=0;j<H->rows();j++)
      if (  H->get(j,i) ) wt++;
    if ( wt < 1)
      cout<<"find column wt < 1 in H. column index "<<i<<" wt = "<<wt<<endl;
  }
  //check row weight of H. check to variable nodes require min wt >=2
  for ( int i =0; i< H->rows();i++){
    int wt=0;
    for ( int j=0;j<H->cols();j++)
      if (  H->get(i,j) ) wt++;
    if ( wt < 2)
      cout<<"find row wt < 2 in H. column index "<<i<<" wt = "<<wt<<endl;
  }
  cout<<H -> rows()<<"  row   rank "<<H->row_rank()<<endl;

  //check min weight of G
  int d = rand_dist(*G);
  cout<< "min weight of G: "<<d<<endl;
  cout<<"finish checking matrices"<<endl;
  return 0;
}

int remove_rows(GF2mat *G, bvec rows_to_remove){
  //remove row in G, corresponding to one in bvec rows_to_remove
  int rows=G->rows(), cols=G->cols();
  ivec perm(rows);
  int t=0;
  //  cout<<perm<<endl;
  for ( int i =0; i<rows;i++){
    if (rows_to_remove(rows-i-1)){
      perm(rows-1-t)= rows-i-1;
      t++;
    }
  }
  //  cout<<"t = "<<t<<endl;
  //cout<<perm<<endl;
  int t1=0;
  for (int i = 0; i< rows;i++){
    if (rows_to_remove(rows-i-1)){
      t1++;
    }
    //perm(rows-i-t-1)=rows-i+t-t1;
    else{
      perm(rows-i-t-1+t1)=rows-i-1;
    }
    // cout<<i<<endl;
    //cout<<perm<<endl;
  }
  G->permute_rows(perm,false);//true or false
  G->set_size(rows-t,cols,true);
  return 0 ;
}

int remove_cols(GF2mat *G, bvec cols_to_remove){
  GF2mat H = G->transpose();
  remove_rows(&H, cols_to_remove);

  GF2mat H1 = H.transpose();
  *G = H1;
  return 0;
}

int remove_cols_mat(mat *G ,bvec cols_to_remove){
  //have to write a seperate one for K
  int size = cols_to_remove.size();
  int t=0;
  for ( int i = size-1; i>-1;i--){
    if ( cols_to_remove(i) ){
      for ( int j=i;j<size-t-1;j++){
	G->swap_cols(j,j+1);
      }
      t++;
    }
  }
  int rows=G->rows(),cols=G->cols();
  G->set_size(rows,cols-t,true);
  return 0;
}

//return true if H has row with weight 1. not in use
bool check_row_weight(GF2mat *H){
  for ( int i =0; i< H->rows();i++){
    int wt=0;
    for ( int j=0;j<H->cols();j++)
      if (  H->get(i,j) ) wt++;
    if ( wt < 2){
      //  cout<<"find row wt < 2 in H. row index "<<i<<" wt = "<<wt<<endl;
      return true;
    }
  }
  return false;
}

//remove rows with weight one in H, and fix other matrices
int reduce_matrices( GF2mat * G, GF2mat * H, GF2mat *U, GF2mat * W, mat * K){
  //now remove weight one rows in H by removing corresponding column in H and other matrices.
  
  //check row weight of H. check to variable nodes require min wt >=2
  //cout<<"start removing columns"<<endl;
  bvec columns_to_remove=zeros_b(H->cols());
  int rt=0;//,ct=0; //rows_to_remove_count=0; columns_to_remove_count=0
  bvec rows_to_remove=zeros_b(H->rows());
  for ( int i =0; i< H->rows();i++){
    int wt=0;
    for ( int j=0;j<H->cols();j++)
      if (  H->get(i,j) ) wt++;
    if ( wt < 2){
      //  cout<<"find row wt < 2 in H. row index "<<i<<" wt = "<<wt<<endl;
      rows_to_remove.set(i,1);
      rt++;
      for ( int j=0;j<H->cols();j++)
	if (  H->get(i,j) ) columns_to_remove.set(j,1);
    }
  }

  if (rt>0){//remove columns if needed
    
    remove_rows(H,rows_to_remove); //try only remove rows
    
    /*
    remove_cols(H,columns_to_remove);
    remove_cols(G,columns_to_remove);
    remove_cols(U,columns_to_remove);
    remove_rows(H,rows_to_remove);
    bvec columns_to_remove_K(1+columns_to_remove.size());
    columns_to_remove_K.set_subvector(1,columns_to_remove);
    columns_to_remove_K.set(0,0);
    remove_cols_mat(K,columns_to_remove_K);
    */
  }
  return 0;
}

//for partial sum decoder
//read those matrices in the folder 
int read_matrices_for_partial_sum(char * filename_prefix, char * filename_suffix, GF2mat * G, GF2mat * H, GF2mat *U, GF2mat * W, mat * K){

  *G = read_one_matrix(filename_prefix,(char *) "g",filename_suffix);
  *H=read_one_matrix(filename_prefix,(char *) "h",filename_suffix);
  *U=read_one_matrix(filename_prefix,(char *) "u",filename_suffix);
  *W=read_one_matrix(filename_prefix,(char *) "w",filename_suffix);
  //read K seperately because it is a mat other than a GF2mat
  char filename[255];
  sprintf( filename,"%s%s%s",filename_prefix, "K", filename_suffix);//append p in the file name
  //  cout<<filename<<endl;
  *K =  MM_to_mat(filename);


  //*H=make_it_full_rank(*H); //not working yet, check it later
  
  //while ( check_row_weight(H) ) { 
  reduce_matrices(G, H,U, W,  K);
  //}
  check_matrices(G, H,U, W,  K);
    
  //  GF2mat G0=MM_to_GF2mat( (char *) "mathematica/matrices/surf5_g_1.mtx" );
  
  return 0;
}

/*
int test_bp_syndrome_llr(GF2mat H,  bvec syndrome, vec & LLRin, vec   & LLRout, int exit_iteration){
  LLRout.set(1,2);
  //  LLRin.set(1,3);
  cout<<LLRin<<endl;
  cout<<LLRout<<endl;
  return 1;
}

bool match_syndrome(GF2mat H, vec LLR, bvec syndrome){
  bvec  error = LLR < 0;
  return GF2mat(H*error - syndrome).is_zero();
}
*/


/*
int bp_syndrome_llr(const GF2mat H,  bvec syndrome,  vec & LLRin, vec   & LLRout, int exit_iteration, int decode_mode){
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

  // *********************************
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
		  str += to_string(llrs(i,k))+",";* /
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
	  }* /	
    }
    if (debug) cout<<"update_count = "<<update_count<<", LLRout = "<<floor(LLRout)<<endl ;
    //if (debug) cout<<"update_count = "<<update_count<<endl;
    //if (debug) draw_toric_x_error(LLRout<0);
    update_count++;
    if ( match_syndrome(H, LLRout, syndrome) ){
      break;
    }        
  }
    
  if (debug) cout<<"LLRout = "<<LLRout<<endl;

  //  if (debug) cout<<"llrs = "<<llrs<<endl;
  //if (debug) cout<<"LLRs = "<<LLRs<<endl;

  //not converge, output negative value
  if (! match_syndrome(H, LLRout, syndrome) ){
    update_count = - update_count;
  }
    
  return update_count ;
}
	*/
