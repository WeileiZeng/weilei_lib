#include "weilei_lib.h"

void test_mmio();

void test_getC();

void test_CSS_code();

void test_product_code();

void test_classical_code();

void test_CSS_Code_IO();

void test_decode();

int main(){
  std::cout<<" --------------------- begin test"<<std::endl;

  //  test_getC();
  //  test_classical_code();
  // test_CSS_code();
  //  test_product_code();
  //  test_mmio();
  //  test_CSS_Code_IO();  //inside test_CSS_code();

  test_decode();
  std::cout<<" --------------------- finish test"<<std::endl;
  return 0;
}





// ============    implementations   ==============



void test_CSS_Code_IO(){
  return;
}

void test_mmio(){
  itpp::GF2mat G = common::get_check_code743(7);
  //  GF2mat_to_MM(G, "tmp/G.mm");
  // itpp::GF2mat H = MM_to_GF2mat("tmp/G.mm");
  //std::cout<<H<<std::endl;
  return;
}

void test_getC(){
  std::cout<<"---------------------- begin test for getC()"<<std::endl;
  ClassicalCode code;
  code.get_743_code(7);
  code.dist();
  itpp::GF2mat H = code.H;
  std::cout<<"H"<<H<<std::endl;

  itpp::GF2mat G = H;//common::nullSpace(H);
  std::cout<<"G"<<G<<std::endl;
  if ( ( G*H.transpose() ).is_zero() ){
    std::cout<<"G*H^T=0\n";
  }
  itpp::GF2mat C = common::getC(G,H);
  std::cout<<"C"<<C<<std::endl;
  return;
}


void test_classical_code(){
  std::cout<<"begin test for ClassicalCode"<<std::endl;
  ClassicalCode code;
 
  if ( code.is_defined ){
    std::cout<<"code is_defined"<<std::endl;
  }else{
    std::cout<<"code is_defined"<<std::endl;
  }

  //Hamming code
  //  code.n = 7;
  code.get_743_code(7);
  code.full_rank();
  code.title="Hamming code";
  code.info();
  code.d = code.dist();
  std::cout<<"The Hamming [7,4,3] code has distance d="<<code.d<<std::endl;

  //repetition code
  code.title="repetition code";
  code.get_repetition_code(7);
  code.full_rank();
  //  code.info();
  code.d = code.dist();
  std::cout<<"repetition code: d="<<code.d<<std::endl;
  std::cout<<code<<std::endl;
  std::cout<<"finish test for ClassicalCode"<<std::endl;
  return;
}


void test_CSS_code(){
  std::cout<<"begin test for CSSCode"<<std::endl;
  CSSCode code;
  code.n = 7;
  code.title="Quantum hamming 713 code";
  code.get_713_code();
  //  std::cout<<code<<std::endl;
  //  code.full_rank();
  //  code.info();
  //  code.d = 
  code.dist();
  std::cout<<"[7,1,3] Hamming code: dx="<<code.dx<<std::endl;
  std::cout<<code<<std::endl;
  std::cout<<"finish test for CSSCode"<<std::endl;
  std::cout<<"begin test for ProductCSSCode"<<std::endl;

  //random CSS code
  CSSCode codeR;
  codeR.title="random code";
  codeR.n=7;
  codeR.Gx_row=3;
  codeR.Gz_row=3;
  codeR.id_Gx=3511;
  codeR.id_Gz=2657;
  codeR.generate_by_id(0);
  codeR.dist();

  std::cout<<codeR<<std::endl;
  std::cout<<codeR.Gx<<std::endl;
  std::cout<<codeR.Gz<<std::endl;


  ProductCSSCode codeP;
  codeP.n=10;
  codeP.title="some Product code codeP";
  std::cout<<codeP<<std::endl;

  ConcatenatedProductCSSCode codeCP;
  codeCP.n=61;
  std::cout<<codeCP<<std::endl;


  SubsystemProductCSSCode codeSP;
  codeSP.n=9;
  SubsystemProductCSSCode codeSP2;
  codeSP=codeSP2;
  std::cout<<codeSP<<std::endl;
  //  codeSP.info();



  std::cout<<"finish test for ProductCSSCode"<<std::endl;


  std::cout<<"start test for CSSCode IO"<<std::endl;

  CSSCode codeS, codeL;
  codeR.save("tmp/testCode");

  codeL.load("tmp/testCode");
  std::cout<<codeL<<std::endl;
  std::cout<<codeL.Gx<<std::endl;
  std::cout<<codeL.Gz<<std::endl;
  std::cout<<"finish test for CSSCode IO"<<std::endl;

  return;
}

void test_product_code(){
  std::cout<<"begin test for product Code"<<std::endl;
  CSSCode code;
  code.n = 7;
  code.title="Quantum hamming 713 code";
  code.get_713_code();
  //  std::cout<<code<<std::endl;
  //  code.full_rank();
  //  code.info();
  //  code.d = 
  code.dist();
  std::cout<<"[7,1,3] Hamming code: dx="<<code.dx<<std::endl;
  std::cout<<code<<std::endl;
  std::cout<<"finish test for CSSCode"<<std::endl;
  std::cout<<"begin test for ProductCSSCode"<<std::endl;

  //random CSS code
  CSSCode codeR;
  codeR.title="random code";
  codeR.n=7;
  codeR.Gx_row=3;
  codeR.Gz_row=3;
  codeR.id_Gx=3511;
  codeR.id_Gz=2657;
  codeR.generate_by_id(0);
  codeR.dist();

  std::cout<<codeR<<std::endl;
  std::cout<<codeR.Gx<<std::endl;
  std::cout<<codeR.Gz<<std::endl;
  //  std::cout<<codeR.Cx<<std::endl;


  CSSCode codeA=code, codeB=codeR;
  codeA.Cx=common::getC(codeA.Gx,codeA.Gz);
  codeA.Cz=common::getC(codeA.Gx,codeA.Gz, 1);
  codeB.Cx=common::getC(codeB.Gx,codeB.Gz);
  codeB.Cz=common::getC(codeB.Gx,codeB.Gz,1);


  ProductCSSCode codeP;
  codeP.n=10;
  codeP.title="some Product code codeP";
  std::cout<<codeP<<std::endl;

  ConcatenatedProductCSSCode codeCP;
  codeCP.n=61;
  std::cout<<codeCP<<std::endl;


  SubsystemProductCSSCode codeSP(codeA,codeB);
  codeSP.title="AB";
  std::cout<<codeSP<<std::endl;
  codeSP.product();

  codeSP.dx = common::quantum_dist_v2((codeSP.Gx), (codeSP.Hz));
  codeSP.dz = common::quantum_dist_v2((codeSP.Gx), (codeSP.Hz));
  std::cout<<"codeSP.dx = "<<codeSP.dx<<std::endl;
  std::cout<<"codeSP.dz = "<<codeSP.dz<<std::endl;
  codeSP.info();

  std::cout<<"finish test for ProductCSSCode"<<std::endl;
  return;
}


//should move to common::
int weight(itpp::bvec &b){//we can use BPSK::count(zero_b(N),b) to count the weight
  int length=b.length();
  int weight=0;
  for (int i=0;i<length;i++){
    if(b.get(i)){
      weight++;
    }
  }
  return weight;
}

//double decode(itpp::bvec e_in, itpp::bvec e_out,double p){

bool decode(itpp::GF2mat Gx, itpp::GF2mat Gz, double p){
  return false;
}

double simulate(itpp::GF2mat Gx, itpp::GF2mat Gz, double p){
  //  itpp::GF2mat Gx,Gz;
  //G and S must be full rank
  itpp::GF2mat G=Gz,S=Gx;

  itpp::RNG_randomize();//get randome seed 
  //#set up parameters#
  int e_try=100;//number of random errors generated
  int perm_try=100;//20;//5;//number of trails of random window / permutation;
  int N=G.cols();///2;//number of qubits, size of the lattice

  itpp::GF2mat T,U;
  itpp::ivec P;//used for all the gaussian elimination in this file, T_fact(T,U,P)

  //cout<<"test commutation: G_t*(S.transpoze()) is "<<G_t*(S.transpose())<<endl;
  //define H=(S_x|S_z|s)
  //find Q, the dual of H, which is e+G+L, HQ^T=0
  
  itpp::GF2mat H=S;    
  H.set_size(S.rows(),S.cols()+1,true);//H=(H_x,s)//H=(H_x,H_z,s)
  G.set_size(G.rows()+1,G.cols(),true);//add one row, use G here because we will use diff=e+X later
 
  itpp::bvec e_t = itpp::zeros_b(N);//e_t(2*N);//e_tilde=e_z//(e_z,e_x)
  //read error from file
  //  GF2mat E_input=MM_to_GF2mat(filename_E);//here the input is the output of the nonconvergent cases after BP decoding. The first row of this vector is an extra zero vector.
  //int e_try=E_input.rows()-1;//number of total input errors
  //   e_try=(10<e_try) ? 10 :e_try; limit it to 10 to smaller
  //GF2mat E_input_good(e_t,false),E_input_bad(e_t,false),E_output_good(e_t,false),E_output_bad(e_t,false);//false for row vectors
  
  itpp::ivec perm(S.cols()+1);
  perm.set(S.cols(),S.cols());//the last col for syndrome, is fixed

  int e_bad=0;//count of bad errors
  for(int i1=0;i1<e_try;i1++){
    for (int i2=0;i2<2*N;i2++){//setup random error with error rate p
      e_t.set(i2,(itpp::randu()-p<0)? 1:0); 
    }

    //    e_t=E_input.get_row(i1+1);//get input error
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
      H.permute_cols(perm, true);//permute back for another run
    }
  
    //check the error and syndrome    
    itpp::bvec diff_t=e_d+e_t;//same format, (e_z,e_x)
      
    // bvec diff=get_tilde(diff_t);
    //      G.set_row(G.rows()-1,diff);//check if the diff belong to gauge group
    G.set_row(G.rows()-1,diff_t);//check if the diff belong to gauge group
    //    std::cout<<G.rows()<<", rank of new G = "<<G.row_rank()<<std::endl;
    if(G.rows()==G.row_rank()){//not belond to gauge -> belongs to logical group -> bad error
      e_bad++;
      //cout<<"BAD* ";
      //cout<<"weight(e_t) = "<<weight(e_t)<<", wmin = "<<wmin<<", weight(diff_t) = "<<weight(diff_t)<<endl;
      //      E_input_bad = append_vector(E_input_bad,e_t);
      //      E_output_bad = append_vector(E_output_bad,e_d);
      
    }else{
      //      std::cout<<e_d<<", "<<e_t<<", "<<diff_t;
      //      std::cout<<"__good__"<<std::endl;
      //      E_input_good = append_vector(E_input_good,e_t);
      //      E_output_good = append_vector(E_output_good,e_d);
      
    }
    //cout<<S*diff_t<<endl;//check if get the same syndrome
  }
  double failure_rate=1.0*e_bad/e_try;

  //cout<<"E_output_bad "<<E_output_bad<<endl;
  //cout<<"p="<<p<<", failure_rate="<<failure_rate<<endl;
  //cout<<"number of input errors:"<<e_try<<endl;
  //cout<<"number of bad errors:"<<e_bad<<endl;

    
  //print result
  std::cout<<"# of bonds/qubits N = "<<N
      <<", # of total input error e_try= "<<e_try
      <<", p = "<<p
      <<", failure_rate = "<<failure_rate
	   <<std::endl;
    
  return failure_rate;
}




void test_decode(){
  std::cout<<"begin test_decode"<<std::endl;
  CSSCode code;
  code.n = 7;
  code.title="Quantum hamming 713 code";
  code.get_713_code();
  //  std::cout<<code<<std::endl;
  //  code.full_rank();
  //  code.info();
  //  code.d = 
  code.dist();
  code.k = code.n - code.Gx.row_rank() - code.Gz.row_rank();
  //std::cout<<"[7,1,3] Hamming code: dx="<<code.dx<<std::endl;
  std::cout<<code<<std::endl;
  std::cout<<"finish generating Steane code"<<std::endl;

  //  double decode(itpp::GF2mat Gx, itpp::GF2mat Gz, double p)
  double p = 0.1;
  code.Gx = common::make_it_full_rank(code.Gx);
  code.Gz = common::make_it_full_rank(code.Gz);
  double p_block = simulate(code.Gx, code.Gz, p);


  std::cout<<"finish test_decode"<<std::endl;
  return;
}
