//Weilei Zeng, April 28
//return min weight of rows in C, which is the distance of the code
//random window method is applied with default number of permutation 10.

//#include "weilei_lib/my_lib.h"
#include "dist.h"
#include "lib.h"
#include "mm_read.h"
#include "mm_write.h"


int min_wt_decoding(itpp::GF2mat C){
  //  std::cout<<C<<endl;
  //when C is not so large, we can apply the true min weight decoding
  //  int rowC=C.rows();
  //  std::cout<<"rowC = "<<rowC<<endl;
  itpp::bvec alpha(C.rows());//a binary vector for linear combination of codeword
  itpp::GF2mat alphaM(1,C.rows());
  int max = pow(2,C.rows());
  itpp::bvec codeword(C.cols());
  int wt, min_wt=C.cols();
  itpp::bvec zero = itpp::zeros_b(C.cols());
  //  std::cout<<"min_wt_decodnig"<<endl;
  //  std::cout<<"max = "<<max<<endl<<C.rows()<<endl;
  //  itpp::bvec min_codeword=zero;
  for ( int i =1;i<max;i++){
    //    itpp::dec2bin(i,alpha);
    alpha = itpp::dec2bin(C.rows(),i);
    alphaM.set_row(0,alpha);
    //std::cout<<i<<std::endl;
    //    std::cout<<"alpha = "<<alpha<<std::endl;
    codeword=(alphaM*C).get_row(0);
    wt = itpp::BERC::count_errors(zero,codeword);
    //    if ( wt < min_wt ) min_codeword=codeword;
    min_wt = (wt<min_wt)? wt : min_wt;
    if ( min_wt == 1) return 1;
  }
  //  std::cout<<"min wt codeword: "<<min_codeword<<std::endl;
  return min_wt;    
}

//G for gauge operators, and C for bare logical operators
//code word c = alpha_C*C+alpha_G*G, where alpha_C \neq 0
int min_wt_decoding(itpp::GF2mat C,itpp::GF2mat G){
  //  std::cout<<"call min_wt_decoding"<<std::endl;
  // make sure G and C are full rank before calling this function. Otherwise it is a waste of computing power.
  int C_rows=C.rows(), G_rows=G.rows();
  int dec_C=(int) pow(2, C_rows);
  int dec_G=(int) pow(2, G_rows);
  int N=C.cols();
  itpp::bvec bvec_C;
  itpp::bvec bvec_zero=itpp::zeros_b(N);
  itpp::GF2mat alpha_C(1,C_rows), alpha_G(1,G_rows);
  int wt=N, min_wt=N;
  //I should save a static copy of this alpha matrix, instead of generate it in two for loops every time.
  for ( int i = 1; i < dec_C ; i++){
    alpha_C.set_row(0,itpp::dec2bin(C_rows,i));
    for ( int j = 0; j < dec_G; j++){
      alpha_G.set_row(0,itpp::dec2bin(G_rows,j));
      bvec_C = (alpha_C * C + alpha_G * G).get_row(0);
      //    std::cout<<bvec_C<<std::endl;
      wt = itpp::BERC::count_errors(bvec_zero,bvec_C);
      //      if ( wt < min_wt ) min_codeword=codeword;
      min_wt = (wt<min_wt)? wt : min_wt;      
      if ( min_wt == 1 ) return 1;
    }
  }
  //  std::cout<<"                          finish min_wt_decoding"<<std::endl;
  return min_wt;
}


int save_dist(int d,char * filename){
  itpp::mat mat_d(1,1);
  mat_d.set(0,0,d);
  mat_to_MM(mat_d,filename);
  return 0;
}



int rand_dist(itpp::GF2mat C, int perm_try){//default perm_try=10
  //  return min wt of rows in C
  if (C.rows()<7){ //for small codes, use true min wt decoding
    return min_wt_decoding(C);
  }
  //use random window decoder to find min wt of C
  //  RNG_randomize(); do not use it here. run it in the main program
  itpp::bvec row_vec,zero=itpp::zeros_b(C.cols());
  int wt,min_wt=C.cols();
  itpp::ivec perm;
  itpp::GF2mat T,U;
  itpp::ivec P;
  for (int j=0;j<perm_try;j++){
    perm = sort_index( itpp::randu(C.cols()) );
    C.permute_cols(perm,false);
    //no need to permute back; can also permute the rows
    //permute rows also
    perm = sort_index( itpp::randu(C.rows()) );
    C.permute_rows(perm,false);
    C.T_fact(T,U,P);
    for (int i = 0;i<C.rows();i++){
      row_vec = C.get_row(i);
      wt = itpp::BERC::count_errors(zero,row_vec);
      if (wt < min_wt){
	min_wt = wt;
	//	std::cout<<row_vec<<std::endl;
      }
    }
  }
  return min_wt;
}
int classical_dist(itpp::GF2mat G){
  //return distance of a classical code GH^t=0
  //G is parity check matrix
  itpp::GF2mat T,U;
  itpp::ivec P;
  int rank_of_G = G.transpose().T_fact(T,U,P);
  if ( rank_of_G == G.cols()  ){
    return INF;//999;//999 for infinity
  }
  itpp::GF2mat H = T.get_submatrix(rank_of_G,0,G.cols()-1,G.cols()-1);
  if (H.rows()<7){//use true min wt decoding for small codes.
    return min_wt_decoding(H);
  }
  int min_wt = rand_dist(H);//default permutation = 10
  return min_wt;
}


//return H such that GH^T = 0, and rank G + rank H = n = full rank
itpp::GF2mat nullSpace(itpp::GF2mat G){
  itpp::GF2mat T,U; itpp::ivec P;
  int n=G.cols();
  int rank_of_G = G.transpose().T_fact(T,U,P);
  //  itpp::GF2matPrint(T,"T");
  itpp::GF2mat Q=T.get_submatrix(rank_of_G,0,n-1,n-1);
  return Q;
}


itpp::GF2mat getC(itpp::GF2mat G_x,itpp::GF2mat G_z,int flip){
  //return C_x
  //flip=1 to get C_z
  if (flip==1){
    itpp::GF2mat temp=G_x;
    G_x=G_z;
    G_z=temp;
  }
  itpp::GF2mat T,U;
  itpp::ivec P;
  int rank_of_G_z =   G_z.transpose().T_fact(T,U,P);
  itpp::GF2mat Q=T.get_submatrix(rank_of_G_z,0,G_z.cols()-1,G_z.cols()-1);//Q include G_x and C_x/L_x
  itpp::GF2mat GQ=G_x.concatenate_vertical(Q);

  GQ.T_fact(T,U,P);
  int rank_of_G_x = G_x.row_rank();
  int rank_of_Q = Q.rows();
  if (rank_of_G_x == rank_of_Q){
    std::cout<<"getC(): It is not a quantum code:zero rank for codeword space"<<std::endl;
    //    return;
  }
  if ( G_x.cols()-rank_of_G_z-rank_of_G_x < 1){
    std::cout<<"empty code space"<<std::endl;
    throw "empty code space";
  }
  //  itpp::GF2matPrint(G_x,"G_x");
  //  itpp::GF2matPrint(U,"U");
  itpp::GF2mat C = U.get_submatrix(rank_of_G_x,0,rank_of_Q-1,G_x.cols()-1 );
  
  C.permute_cols(P,true);//codewords/logical group 
  //check  if ((G_z*C.transpose()).is_zero() ){    std::cout<<"GOOD C"<<std::endl;  }
  //  std::cout<<"get C"<<std::endl;
  return C;
}


int quantum_dist_v2(itpp::GF2mat G_x, itpp::GF2mat G_z, int flip){//without expected value
  //right or x  distance of (G_x,G_z)
  //flip left and right if flip = 1;
  int trialQ=500;//100 is good for not so big codes;permute GQ this max amount of time
  //  int trialQflag=1;//a falg to control the max amout of permutation
  
  if (flip==1){//flip G_x and G_z
    itpp::GF2mat temp=G_x;    G_x=G_z;    G_z=temp;
  }



  itpp::GF2mat T,U;  itpp::ivec P;
  int rank_of_G_z =   G_z.transpose().T_fact(T,U,P);

  // for small code, use min_wt_decoding to return X distance
  if (G_z.cols() - rank_of_G_z < 11){
    return min_wt_decoding(getC(G_x, G_z), G_x);
  } 


  itpp::GF2mat Q=T.get_submatrix(rank_of_G_z,0,G_z.cols()-1,G_z.cols()-1);//Q include G_x and C_x/L_x
  itpp::GF2mat GQ=G_x.concatenate_vertical(Q);
  int min_wt=GQ.cols(),wt;

  for ( int iq=0;iq<trialQ;iq++){
    itpp::ivec perm = sort_index(  itpp::randu( GQ.cols()  ));//random permutation
    GQ.permute_cols(perm,false);
    GQ.T_fact(T,U,P);
    int rank_of_G_x = G_x.row_rank();
    int rank_of_Q = Q.rows();

    if (rank_of_G_x == rank_of_Q){
      return INF;//999 for infinity
    }
    itpp::GF2mat C = U.get_submatrix(rank_of_G_x,0,rank_of_Q-1,G_x.cols()-1 );
    C.permute_cols(P,true);//codewords/logical group //not necessary to permute it back here
  
    wt = rand_dist(C);//defauylt permutation = 10
    trialQ=(wt<min_wt)? ( (10*iq > trialQ)? 10*iq : trialQ ):trialQ;//make sure this is the true min weight
    min_wt=(wt<min_wt)? wt:min_wt;
    if (min_wt ==1) return 1;

    //    std::cout<<"iq = "<<iq<<", [wt="<<wt<<"] "<<std::endl;;
    //  std::cout<<"got min wt of logical operator C  = "<<min_wt<<std::endl;
    //save_dist(min_wt,filename_dist);
  }
  return min_wt;
}
int quantum_dist(itpp::GF2mat G_x, itpp::GF2mat G_z, int dist_expected, int debug, int flip){
  //right or x  distance of (G_x,G_z)
  //flip left and right if flip = 1;
  int trialQ=50000;//1000;permute GQ this max amount of time
  if ( dist_expected > 10 ) trialQ = trialQ*2;
  int trialQflag=1;//a flag to adjust the max amount of permutation
  
  if (flip==1){//flip G_x and G_z
    itpp::GF2mat temp=G_x;    G_x=G_z;    G_z=temp;
  }

  itpp::GF2mat T,U;  itpp::ivec P;
  int rank_of_G_z =   G_z.transpose().T_fact(T,U,P);
  itpp::GF2mat Q=T.get_submatrix(rank_of_G_z,0,G_z.cols()-1,G_z.cols()-1);//Q include G_x and C_x/L_x

  //  std::cout<<Q.row_rank()<<","<<G_z.row_rank()<<","<<G_z.cols()<<std::endl;

  itpp::GF2mat GQ=G_x.concatenate_vertical(Q);
  int min_wt=GQ.cols(),wt;

  for ( int iq=0;iq<trialQ;iq++){
    itpp::ivec perm = sort_index(  itpp::randu( GQ.cols()  ));//random permutation
    GQ.permute_cols(perm,false);
    GQ.T_fact(T,U,P);
    int rank_of_G_x = G_x.row_rank();
    int rank_of_Q = Q.rows();

    if (rank_of_G_x == rank_of_Q){
      return INF;//999 for infinity
    }
    itpp::GF2mat C = U.get_submatrix(rank_of_G_x,0,rank_of_Q-1,G_x.cols()-1 );
    C.permute_cols(P,true);//codewords/logical group //not necessary to permute it back here
    //Question 1: does the row in C include some stabilizer generators which may increase its weight?
  
    wt = rand_dist(C);//defauylt permutation = 10
    min_wt=(wt<min_wt)? wt:min_wt;
    //std::cout<<"iq = "<<iq<<", [wt="<<wt<<"] "<<std::endl;;
    //  std::cout<<"got min wt of logical operator C  = "<<min_wt<<std::endl;
    //save_dist(min_wt,filename_dist);
    int max_trial = 0;//no need to run C again, always get the same result
    for (int i =0;i<max_trial;i++){
	  wt=rand_dist(C);
	  min_wt=(wt<min_wt)? wt:min_wt;
    }
    if (trialQflag) {//adjust the max number of permutation, only do this once
        if (min_wt <= dist_expected){
	  trialQflag=0;
	  trialQ = 10*iq;
	  // continue to run to see if smaller distance can be achieved.
	  trialQ=(trialQ<1000)? 1000:trialQ;
	  if (debug) std::cout<<"quantum_dist: reach min distance when iq = "<<iq<<", continue to run with trialQ = "<<trialQ<<std::endl;
	}
    }
  }
  return min_wt;
}

int hypergraph_dist(itpp::GF2mat Aj, itpp::GF2mat Ajplus,int dist_expected,int flip){
  //left distance of (Aj,Ajplus^T)
  //right distance of (G_x,G_z)
  //flip left and right if flip = 1
  
  itpp::GF2mat G_z = Aj;
  itpp::GF2mat G_x = Ajplus.transpose();
  //check commutation
  /*  if ( (Aj*Ajplus).is_zero()){
    //pass
  }else{
    std::cout<<"It is not a CSS code!"<<std::endl;
    }*/
  
  if (flip==1){
    itpp::GF2mat temp=G_x;
    G_x=G_z;
    G_z=temp;
  }
  
  /* char * filename_G_x = argv[1];
  char * filename_G_z = argv[2];
  char * filename_C_x = argv[3];
  char * filename_dist = argv[4];
  */
  //  itpp::GF2mat G_z = MM_to_GF2mat(filename_G_z);
  itpp::GF2mat T,U;
  itpp::ivec P;
  int rank_of_G_z =   G_z.transpose().T_fact(T,U,P);
  itpp::GF2mat Q=T.get_submatrix(rank_of_G_z,0,G_z.cols()-1,G_z.cols()-1);//Q include G_x and C_x/L_x

  //  itpp::GF2mat G_x=MM_to_GF2mat(filename_G_x);


  itpp::GF2mat GQ=G_x.concatenate_vertical(Q);
  int min_wt=GQ.cols(),wt;
  int trialQ=1000;//1000;
  for ( int iq=0;iq<trialQ;iq++){
    itpp::ivec perm = sort_index(  itpp::randu( GQ.cols()  ));
    GQ.permute_cols(perm,false);
    GQ.T_fact(T,U,P);
    int rank_of_G_x = G_x.row_rank();
    int rank_of_Q = Q.rows();

    if (rank_of_G_x == rank_of_Q){
      //      std::cout<<"It is not a quantum code!"<<std::endl;
      if (dist_expected != 999){
	std::cout<<std::endl<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!NOTICE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!NOTICE!!!!!!!!!!!!!"<<std::endl;
	std::cout<<"("<<dist_expected<<")";
      }
      return 999;//999 for infinity
    }
    itpp::GF2mat C = U.get_submatrix(rank_of_G_x,0,rank_of_Q-1,G_x.cols()-1 );
    C.permute_cols(P,true);//codewords/logical group //not necessary to permute it back here

    //Question 1: does the row in C include some stabilizer generators which may increase its weight?
  
    //itpp::GF2mat_to_MM(C,filename_C_x);
    wt = rand_dist(C);//defauylt permutation = 10
    min_wt=(wt<min_wt)? wt:min_wt;
    //  std::cout<<"got min wt of logical operator C  = "<<min_wt<<std::endl;
    //save_dist(min_wt,filename_dist);
    //    std::cout<<"<"<<wt<<">";      
    int max_trial = 0;//no need to run C again, always get the same result
    for (int i =0;i<max_trial;i++){
        if (min_wt == dist_expected){
	  //great job! return the result
	  break;
	}else if (min_wt < dist_expected){
	  //      std::cout<<"Damn! how could this happen!"<<std::endl;
	  //    std::cout<<"("<<dist_expected<<")";
	  break;
	}else{
	  //continue another run
	  //      std::cout<<"#RUN AGAIN#"<<min_wt;
	  wt=rand_dist(C);
	  min_wt=(wt<min_wt)? wt:min_wt;
	  //	  std::cout<<"<"<<wt<<">";
	}
    }

    if (min_wt == dist_expected){
      break;
    }
  }

  //final check
  if (min_wt != dist_expected){
    std::cout<<std::endl<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!NOTICE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!NOTICE!!!!!!!!!!!!!"<<std::endl;
    std::cout<<"("<<dist_expected<<")";
    
  }


  return min_wt;
}




int draw_toric_x_error(itpp::bvec error_bits){
  //  std::cout<<"draw with weight = "<<weight(error_bits)<<std::endl;
  //  if (true) return 0 ;
  
  //draw ( but not print) error for the toric code
  //in stabilizer generating matrix, we have vertex for X and plaquette for Z
  //use plaquette to check X error and vertex to check Z errors
  //for indexing, start from the plaqutte in the top left. Its top bond has index one, and then move to the right and then move down. After finishing all horizontal bonds in the lattice, Start from its left bond with index n*n, the move right and move down.
  int N=error_bits.length(); //size of the code N=2*n*n
  int n=(int) sqrt(N/2); //size of the lattice

  //first row of horizontal bonds
  std::cout<<" ";
  for ( int j=0; j<n; j++){
    if (error_bits.get(n*0+j)){
      std::cout<<color_text("_ ");
    }else{
      std::cout<<"_ ";
    }
  }
  std::cout<<std::endl;
  for (int i=1;i<n;i++){
    for ( int j=0; j<n; j++){
      if (error_bits.get(n*(i-1)+n*n+j)){
	std::cout<<color_text("|");//"1";//for error bits
      }else{
	std::cout<<"|";
      }
      if (error_bits.get(n*i+j)){
	std::cout<<color_text("_");
	//for error bits
      }else{
	std::cout<<"_";
      }
    }
    std::cout<<std::endl;
  }
  //the last row of vertical bonds
  for ( int j=0; j<n; j++){
    if (error_bits.get(n*(n-1)+n*n+j)){
      //      std::std::cout<<"\033[0;32m 1 \033[0m";
      std::cout<<color_text("| ");
    }else{
      std::cout<<"| ";
    }
  }
  std::cout<<std::endl;
  return 0;
}

int draw_toric_x_error(itpp::bvec error_bits, std::string header){
  std::cout<<header.c_str()<<std::endl;
  draw_toric_x_error(error_bits);
  return 0;
  
}

itpp::bvec find_error(itpp::bvec e_in, itpp::GF2mat H){
  //input: original error and parity check matrix
  //output: an error with same syndrome
  //for principle, see random window decoder


  //  std::cout<<H.rows()<<"  row rank "<< H.row_rank()<<std::endl;
  H = make_it_full_rank(H);  //H may not have full rank
  itpp::GF2mat H0=H;
  itpp::bvec e_t = e_in;//change name
  itpp::bvec s=H*e_t;//syndrome

  H.set_size(H.rows(),H.cols()+1,true);//H=(H_x,s)
  H.set_col(H.cols()-1,s);//add syndrome

  itpp::GF2mat T,U;
  itpp::ivec P;
  H.transpose().T_fact(T,U,P);
  itpp::GF2mat Q=T.get_submatrix(H.rows(),0,H.cols()-1,H.cols()-1);
  itpp::bvec X_t=Q.get_row(Q.rows()-1);//the error detected, X_t=(X_z,X_x,1)
  X_t.set_size(X_t.size()-1,true);
  return X_t;
}


itpp::GF2mat get_check_code734(int L){//L=7n
  //return check matrix code code [7,3,4], find definition in research note.pdf
  itpp::GF2mat G(L,L);
  for ( int i=0;i<L;i++){
    G.set(i,i,1);
    if (i+2>L-1) {
      G.set(i,i+2-L,1);
    } else {
      G.set(i,i+2,1);
    }
    if (i+3>L-1) {
      G.set(i,i+3-L,1);
    }else{
      G.set(i,i+3,1);
    }
  }
  return G;
}

itpp::GF2mat get_check_code743(int L){//L=7n
  //return check matrix code code [7,4,3], find definition in research note.pdf
  itpp::GF2mat G(L,L);
  for ( int i=0;i<L;i++){
    G.set(i,i,1);
    if (i+2>L-1) {
      G.set(i,i+2-L,1);
    } else {
      G.set(i,i+2,1);
    }
    if (i+3>L-1) {
      G.set(i,i+3-L,1);
    }else{
      G.set(i,i+3,1);
    }
    if (i+4>L-1) {
      G.set(i,i+4-L,1);
    }else{
      G.set(i,i+4,1);
    }
  }
  return G;
}
itpp::GF2mat get_check_rept(int L){//return circulant check matrix for repetition code of length L
  itpp::GF2mat a(L,L);
  for (int i=0;i<L-1;i++){//construct a : circulant check matrix for repetition code
    a.set(i,i,1);
    a.set(i,i+1,1);
  }
  a.set(L-1,L-1,1);
  a.set(L-1,0,1);
  //std::cout<<"circulant check matrix for repetition code : a = "<<a<<std::endl;
  return a;
}
itpp::GF2mat get_check(int generator_flag, int L){
  //return check matric a for generating toric code, cubic code and hypercubic code.
  switch(generator_flag){
  case 1: return get_check_rept(L);break;
  case 2: return get_check_code734(L);break;
  case 3: return get_check_code743(L);break;
  }
  //default
  return get_check_rept(L);
  //  return get_check_734(L);//code [7,3,4]
  //  return get_check_743(L);//code [7,4,3]
  //  return get_check_rept(L); //circulant check matrix for repetition code.
}
/*
itpp::LDPC_Code get_test_LDPC_Code(){
  //convert itpp::GF2mat saved in .mm file to LDPC_Code
  //  itpp::GF2mat G=MM_to_GF2mat(filename);
  itpp::GF2mat G = get_check_code734(7*2);
  GF2mat_sparse Gs=G.sparsify();
  itpp::GF2mat_sparse_alist Gsa;
  Gsa.from_sparse(Gs);
  itpp::LDPC_Parity H(Gsa);
  itpp::LDPC_Code C(&H);
  return C;
  }*/

itpp::LDPC_Code GF2mat_to_LDPC_Code(itpp::GF2mat G){
  itpp::GF2mat_sparse Gs=G.sparsify();
  itpp::GF2mat_sparse_alist Gsa;
  Gsa.from_sparse(Gs);
  itpp::LDPC_Parity H(Gsa);
  itpp::LDPC_Code C(&H);
  return C;
}
  
itpp::LDPC_Code MM_to_LDPC_Code(char * filename){
  //convert itpp::GF2mat saved in .mm file to itpp::LDPC_Code
  itpp::GF2mat G=MM_to_GF2mat(filename);
  return  GF2mat_to_LDPC_Code(G);
  /*
  
  itpp::GF2mat_sparse Gs=G.sparsify();
  itpp::GF2mat_sparse_alist Gsa;
  Gsa.from_sparse(Gs);
  itpp::LDPC_Parity H(Gsa);
  itpp::LDPC_Code C(&H);
  return C;*/
}

