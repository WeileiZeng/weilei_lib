//#ifndef CONCATENATION_LIB_H
//#define CONCATENATION_LIB_H
#ifndef PRODUCT_LIB_H
#define PRODUCT_LIB_H

#include "dist.h"
#include <itpp/itbase.h>
//#include <itpp/itcomm.h>
//#include <stdio.h>
#include "weilei_lib.h"


const int MAX_M=6;//maximum of the length of the complex chain
//const int INF=999;//infinity distance


/** \class ClassicalCode
 * a classical binary code
 */
class ClassicalCode{
public:
  itpp::GF2mat G; ///< codeword generating matrix
  itpp::GF2mat H; ///< parity check matrix
  int n; ///< Number of bits
  int k ///< number of encoded bits
  int d; ///< distance 
  int is_defined=0; ///< if G and H has been defined

  //constructor
  ClassicalCode();
  ClassicalCode(itpp::GF2mat G, itpp::GF2mat H);

  //distance estimator
  int dist();
  int min_weight_dist();
  int rand_dist();


  //function
  void info();
  /** return dual code */
  ClassicalCode dual();
  void full_rank();

  //generate sample code
  void get_repetition_code();
  void get_743_code();
  void get_734_code();
};



/** a wrapper of data for a CSS code */
class CSSCode{
public:
  itpp::GF2mat Gx; /* X type parity check matrix */
  itpp::GF2mat Gz;
  itpp::GF2mat Cx;
  itpp::GF2mat Cz;
  itpp::bvec min_weight_codeword;
  int n;
  int Gx_row, Gz_row;
  int id_Gx, id_Gz;   /** id used when enumerating all cases*/

  int dx,dz;
  int is_defined=0, is_C_defined=0;
  
  //constructor
  CSSCode();
  /**
   *@param id_Gax see definition in generate_code()
   *@param id_Gaz see definition in generate_code()
   */
  CSSCode(int na, int Gax_row, int id_Gax, int Gaz_row, int id_Gaz);

  //generating code
  int generate_by_id(int debug);
  int getRandomCode();
  int getGoodCode(int debug);

  //sanity check
  bool is_valid();

  //distance estimation
  /** call rand_dist to estimate dx and dz*/
  void dist();
  int min_weight_dist_x();
  int min_weight_dist_z();
  int rand_dist_x();
  int rand_dist_z();

  //decoding

  //generate sample code
  /** generate 7 qubit hamming code*/
  void get_713_code();
};
 
class ProductCSSCode: public CSSCode{
public:
  CSSCode codeA, codeB;
  ProductCSSCode(){}
  ProductCSSCode(CSSCode codeA_temp, CSSCode codeB_temp);
};


/** a wrapper of data for a product of two CSS codes. */
class SubsystemProductCSSCode : public ProductCSSCode {
public:
  //  string title_str, string note, int mode, int sub_mode_A, int sub_mode_B,     //general info
  // int n_low, int n_high, int k_low, int k_high, int debug,                     //for random simulation
  //  int na;
  //  int Gax_row; int id_Gax; int Gaz_row; int id_Gaz;   //for enumarating all cases
  //  int Gbx_row; int id_Gbx; int Gbz_row; int id_Gbz;   //for enumarating all cases
  //  int is_defined=0;

  //  CSSCode codeA, codeB;
  //  SubsystemProductCode();
  SubsystemProductCSSCode(){}
  SubsystemProductCSSCode(CSSCode codeA_temp, CSSCode codeB_temp){
    codeA=codeA_temp;
    codeB=codeB_temp;
    //    std::cout<<" get codeA with codeA.n = "<< codeA_temp.n<<std::endl;
    if ( codeA.is_defined && codeB.is_defined ){
      //      std::cout<<"both code A and code B are defined"<<std::endl;
      is_defined=1;
    }
  }
};





int getRandomQuantumCode(int n,int Gx_row,int Gz_row, itpp::GF2mat &Gx,itpp::GF2mat &Gz, itpp::GF2mat &Cx,itpp::GF2mat &Cz);

int getGoodQuantumCode(int n,int Gx_row,int Gz_row, itpp::GF2mat &Gx,itpp::GF2mat &Gz, itpp::GF2mat &Cx,itpp::GF2mat &Cz, int debug);

void set_submatrix(itpp::GF2mat & G, itpp::GF2mat sub, int row, int col);

/** A CSS code can be uniquely defined by its dimension and ID
 */
int generate_code(itpp::GF2mat & Gax, itpp::GF2mat & Gaz, int na, int Gax_row, int id_Gax, int Gaz_row, int id_Gaz, int debug);
//moved here cause it is used in CSSCode


/** a wraper */
int generate_code(CSSCode & code, int debug);
//  return generate_code(code.Gx, code.Gz, code.n, code.Gx_row, code.id_Gx, code.Gz_row, code.id_Gz, debug);
//*/





//files used in concatenated codes and product codes.

//int reduce(itpp::GF2mat Gax, itpp::GF2mat Gaz, itpp::GF2mat Gbx, itpp::GF2mat Gbz,int ddax,int ddaz,int ddbx,int ddbz);
//int concatenate(itpp::GF2mat Gax, itpp::GF2mat Gaz, itpp::GF2mat Gbx, itpp::GF2mat Gbz,int ddax,int ddaz,int ddbx,int ddbz);



// a version include both reduce and concatenation
// mode=1 for reduce/subsystem product
// mode=2 for concatenation
//only dz is checked cause dx is known to be tight
int product(itpp::GF2mat Gax, itpp::GF2mat Gaz, itpp::GF2mat Gbx, itpp::GF2mat Gbz,int ddax,int ddaz,int ddbx,int ddbz, int debug, int mode);




#endif //PRODUCT_LIB_H
