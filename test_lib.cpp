#include "weilei_lib.h"

void test_mmio();

void test_getC();

void test_CSS_code();

void test_product_code();

void test_classical_code();

void test_CSS_Code_IO();

int main(){
  std::cout<<" --------------------- begin test"<<std::endl;

  //  test_getC();
  //  test_classical_code();
  // test_CSS_code();
  test_product_code();
  //  test_mmio();
  //  test_CSS_Code_IO();  //inside test_CSS_code();
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
  //  codeSP.n=codeSP.Gx.cols();
  //  itpp::GF2mat Gx = common::make_it_full_rank(codeSP.Gx);
  //  std::cout<<codeSP.Gx.rows()<<Gx.rows()<<std::endl;
  //  std::cout<<codeSP.is_valid()<<std::endl;

  codeSP.dx = common::quantum_dist_v2((codeSP.Gx), (codeSP.Hz));
  codeSP.dz = common::quantum_dist_v2((codeSP.Gx), (codeSP.Hz));
  std::cout<<"codeSP.dx = "<<codeSP.dx<<std::endl;
  std::cout<<"codeSP.dz = "<<codeSP.dz<<std::endl;
  codeSP.info();



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
