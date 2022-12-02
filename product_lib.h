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
 * a classical binary code. It is similar to itpp::LDPC_Code
 * Note not all function are robust against dimension and rank. Those are not implemented for faster speed. One should do sanity check if needed.
 */
class ClassicalCode{
public:
  itpp::GF2mat G; ///< codeword generating matrix
  itpp::GF2mat H; ///< parity check matrix
  int n=-1; ///< Number of bits
  int k=-1; ///< number of encoded bits
  int d=-1; ///< distance 
  int is_defined=0; ///< if G and H has been defined
  std::string title="no title";
  std::string type="ClassicalCode";
  //constructor
  ClassicalCode();
  /**
   *@param G codeword generating matrix
   *@param H parity check matrix
   */
  ClassicalCode(itpp::GF2mat G, itpp::GF2mat H);

  //distance estimator
  /**@returns distance of the code */
  int dist();
  /** min weight decoder */
  int min_weight_dist();
  /** random window decoder */
  int rand_dist();

  //function
  /** print basic infomation */
  void info();
  friend std::ostream& operator<<(std::ostream& os, const ClassicalCode& code);
  /** return dual code */
  ClassicalCode dual();
  /** make G and H full rank */
  void full_rank();

  //generate sample code
  /** generate repetition code*/
  void get_repetition_code(int L);
  /** generate [7,4,3] code according to code.n, which needs to be a multiple of 7*/
  void get_743_code(int L);
  /** generate [7,3,4] code according to code.n, which needs to be a multiple of 7*/
  void get_734_code(int L);

};



/** \class CSSCode
 *a wrapper of data for a CSS code */
class CSSCode{
public:
  itpp::GF2mat Gx; ///< X type parity check matrix 
  itpp::GF2mat Gz; ///< Z type parity check matrix 
  itpp::GF2mat Cx; ///< X type codeword generating matrix
  itpp::GF2mat Cz; ///< Z type codeword generating matrix
  itpp::bvec min_weight_codeword_x;  ///< If decoded, the min weight codeword of type X will be saved 
  itpp::bvec min_weight_codeword_z;  ///< If decoded, the min weight codeword of type Z will be saved 
  int n=-1; ///< code size, number of qubits
  int k=-1; ///< encoded qubits
  int Gx_row, Gz_row;
  int id_Gx;   ///< id used when enumerating all cases, see definition in generate_by_id()
  int id_Gz;   ///< id used when enumerating all cases, see definition in generate_by_id()
  std::string title="no title";
  std::string type="CSSCode";
  int d=-1;   ///< d=min(dx,dz)
  int dx=-1;  ///< X type min distance
  int dz=-1;  ///< Z type min distance
  int is_defined=0, is_C_defined=0;
  
  //constructor
  CSSCode();
  /**
   *@param id_Gax see definition in generate_code()
   *@param id_Gaz see definition in generate_code()
   */
  CSSCode(int na, int Gax_row, int id_Gax, int Gaz_row, int id_Gaz);

  //generate code
  /**generating code by unique id*/
  int generate_by_id(int debug);
  /** get a random code*/
  int getRandomCode();
  /** get a random code with distance > 1*/
  int getGoodCode(int debug);
  /** calculate Cx and Cz based on Gx and Gz */
  int set_up_CxCz();

  //sanity check
  bool is_valid();   ///< check if the code is valid
  void full_rank();  ///< make matrices full rank. Not implemented yet
  void info();   ///< not in use. similar to operator<<()
  friend std::ostream& operator<<(std::ostream& os, const CSSCode& code);  ///< print basic info


  //distance estimation
  /** call rand_dist to estimate dx and dz*/
  void dist();
  int min_weight_dist_x();
  int min_weight_dist_z();
  int rand_dist_x();
  int rand_dist_z();

  //decoding
  /** decode error e 
      
   */
  //  void decode(itpp::bvec e_in, itpp::bvec e_out);
  bool decode(itpp::bvec e_t, const int perm_try, const int debug=0);
  //should add e_out here
  //  double simulate(double p, const int e_try = 1000);
  double simulate(double p, const int e_try = 1000, const int num_cores=16, const int debug = 0);

  //generate sample code
  /** generate 7 qubit hamming code*/
  void get_713_code();

  //I/O
  
  int save(std::string filename_prefix); ///<save code matrices into MM 
  int load(std::string filename_prefix);
};
 
class ProductCSSCode: public CSSCode{
public:
  CSSCode codeA, codeB;

  std::string type="ProductCSSCode";
  ProductCSSCode(){
  }
  ProductCSSCode(CSSCode codeA_temp, CSSCode codeB_temp);
  //shourld I use  virtual  here or not?
  void product(); ///< generate the product code from codeA and codeB, to be impelmented in each derived class

  friend std::ostream& operator<<(std::ostream& os, const ProductCSSCode& code);
};


/** a wrapper of data for a product of two CSS codes. */
class SubsystemProductCSSCode : public ProductCSSCode {
public:
  //  string title_str, string note, int mode, int sub_mode_A, int sub_mode_B,     //general info
  // int n_low, int n_high, int k_low, int k_high, int debug,                     //for random simulation
  itpp::GF2mat Hx;
  itpp::GF2mat Hz;
  std::string      type="SubsytemProductCSSCode";
  //  SubsystemProductCode();
  SubsystemProductCSSCode(){}
  SubsystemProductCSSCode(CSSCode codeA_temp, CSSCode codeB_temp):ProductCSSCode( codeA_temp, codeB_temp){}
  void product();

  friend std::ostream& operator<<(std::ostream& os, const SubsystemProductCSSCode& code);
};

class ConcatenatedProductCSSCode: public ProductCSSCode {
public: 
  std::string     type="ConcatenatedProductCSSCode";
  ConcatenatedProductCSSCode(){
}
  ConcatenatedProductCSSCode(CSSCode codeA_temp, CSSCode codeB_temp):ProductCSSCode( codeA_temp, codeB_temp){}
  //  void product();

  friend std::ostream& operator<<(std::ostream& os, const ConcatenatedProductCSSCode& code);
};


//   ------------------------ finish class definition ----------------------


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



/** a version include both reduce and concatenation
 *@param mode=1 for reduce/subsystem product
 *@param mode=2 for concatenation
 *only dz is checked cause dx is known to be tight
*/
int product(itpp::GF2mat Gax, itpp::GF2mat Gaz, itpp::GF2mat Gbx, itpp::GF2mat Gbz,int ddax,int ddaz,int ddbx,int ddbz, int debug, int mode);




#endif //PRODUCT_LIB_H
