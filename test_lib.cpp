#include "weilei_lib.h"


void test_classical_code(){
  std::cout<<"begin test for ClassicalCode"<<std::endl;
  ClassicalCode code;
 
  if ( code.is_defined ){
    std::cout<<"code is_defined"<<std::endl;
  }else{
    std::cout<<"code is_defined"<<std::endl;
  }

  code.n = 7;
  code.get_743_code();
  code.full_rank();
  code.info();
  code.d = code.dist();
  std::cout<<"[7,4,3] code: d="<<code.d<<std::endl;

  code.get_repetition_code();
  code.full_rank();
  code.info();
  code.d = code.dist();
  std::cout<<"repetition code: d="<<code.d<<std::endl;
  std::cout<<"finish test for ClassicalCode"<<std::endl;
  return;
}


void test_CSS_code(){
  std::cout<<"begin test for CSSCode"<<std::endl;
  CSSCode code;
 
  code.n = 7;
  code.get_713_code();
  //  code.full_rank();
  //  code.info();
  //  code.d = 
  code.dist();
  std::cout<<"[7,1,3] Hamming code: dx="<<code.dx<<std::endl;

  std::cout<<"finish test for CSSCode"<<std::endl;
  return;
}



int main(){
  std::cout<<"begin test"<<std::endl;

  test_classical_code();
  test_CSS_code();

  std::cout<<"finish test"<<std::endl;
  return 0;
}
