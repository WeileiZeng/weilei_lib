#include "weilei_lib.h"


void test_classical_code(){
  std::cout<<"begin test for ClassicalCode"<<std::endl;
  ClassicalCode code;
 
  if ( code.is_defined ){
    std::cout<<"code is_defined"<<std::endl;
  }else{
    std::cout<<"code is_defined"<<std::endl;
  }

  //Hamming code
  code.n = 7;
  code.get_743_code();
  code.full_rank();
  code.title="Hamming code";
  code.info();
  code.d = code.dist();
  std::cout<<"The Hamming [7,4,3] code has distance d="<<code.d<<std::endl;

  //repetition code
  code.title="repetition code";
  code.get_repetition_code();
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

  ProductCSSCode codeP;
  codeP.n=10;
  codeP.title="some codeP";
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
  return;
}



int main(){
  std::cout<<"begin test"<<std::endl;

  test_classical_code();
  test_CSS_code();

  std::cout<<"finish test"<<std::endl;
  return 0;
}
