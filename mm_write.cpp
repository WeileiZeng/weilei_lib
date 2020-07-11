/** \file mm_write.cpp
 *\author Weilei Zeng
 *   revised from Matrix Market I/O example program
 *   (See http://math.nist.gov/MatrixMarket for details.)
 */
//recomend file format for Matrix Market file is .mm or .mtx
//the segmentation fault is finally solved by reserve memory for I,J and Vval. May 2018
//trouble shooting: segmentation fault:  directory not exist; wrong folder name

#include "mmio.h"
//#include <stdio.h>
//#include <stdlib.h>
//#include <fstream> //ofstream ifstream
#include <itpp/itbase.h>


//int GF2mat_to_MM(GF2mat G, char* file_name="mm_temp.dat")
int GF2mat_to_MM(itpp::GF2mat G, char* file_name, int debug)
{
  int M=G.rows(),N=G.cols(),nt=M*N;//nt is the total number of elements, will find nz later. nz is the number of non-zero elements
  //make nt smaller; the size of int[] should be less than about 2000000+; otherwise it creat segmentation fault
  //this could be fixed by seperate int[] into smaller ones, but in this case the matrix is already too bigfor any calculation
  // std::cout<<"density of GF2mat: "<<G.density()<<endl;
  nt=floor(1+nt*(G.density()+0.001)  );

  int *I, *J;
  I = (int *) malloc(nt * sizeof(int));
  J = (int *) malloc(nt * sizeof(int));
  double * val;
  val = (double *) malloc(nt * sizeof(double));
  
  //  int I[nt],J[nt];//the location of nonzero elements
  //usually nt is much larger than nz, where nz/nt is the density of the code. There is some wasted memory here.
  //  double val[nt];//have to use double to save the binary value in order to match the format of the file, didn't affect the result so far
  int pos=0;
  for(int m=0;m<M;m++){
    //      std::cout<<"\t m="<<m<<",";
      for(int n=0;n<N;n++){
	  if(G.get(m,n)){//if it is nonzero
	    I[pos]=m;
	    J[pos]=n;
	    val[pos]=1;//we know it is one. so dont use G.get(m,n)
	    pos++;
	  }
      }
    }
    int nz=pos;//number of nonzero elements
    //std::cout<<"number of ones in the matrix is nz= "<<nz<<std::endl;
    MM_typecode matcode;
    mm_initialize_typecode(&matcode);
    mm_set_matrix(&matcode);
    mm_set_coordinate(&matcode);
    mm_set_real(&matcode);

    //create a file and wrote into it.
    FILE *fout;
    fout=fopen(file_name,"w");
    mm_write_banner(fout, matcode); 
    mm_write_mtx_crd_size(fout, M, N, nz);

    /* NOTE: matrix market files use 1-based indices, i.e. first element
      of a vector has index 1, not 0.  */

    //it is not necessary to seperate the reading and writing in 2 loops
    for (int i=0; i<nz; i++)
        fprintf(fout, "%d %d %10.3g\n", I[i]+1, J[i]+1, val[i]);

    fclose(fout);
    if ( debug )   std::cout<<"wrote the matrix (density:"<<G.density()<<") into file "<<file_name<<std::endl;
    return 0;
	
}

//int mat_to_MM(mat G, char* filename="mm_temp.dat")
int mat_to_MM(itpp::mat G, char * filename){//use this when G is not a binary matrix
  //dose not work sometimes and not sure why. Segmentation fault in line:   mm_write_banner(fout, matcode);
  int M=G.rows(),N=G.cols(),nt=M*N;//nt is the total number of elements, will find nz later. nz is the number of non-zero elements
  //make nt smaller; the size of int[] should be less than about 2000000+
  //nt=floor(1+nt*G.density());//design to be not sparse and small size
  MM_typecode matcode;
  int I[nt],J[nt];//the location of nonzero elements
  //usually nt is much larger than nz, where nz/nt is the density of the code. There is some wasted memory here.
  double val[nt];//have to use double to save the binary value in order to match the format of the file, didn't affect the result so far
  int pos=0;
  for(int m=0;m<M;m++){
    for(int n=0;n<N;n++){
      if(G.get(m,n)){//if it is nonzero
	I[pos]=m;
	J[pos]=n;
	val[pos]=G.get(m,n);
	pos++;
      }
    }
  }
  int nz=pos;//number of nonzero elements
  mm_initialize_typecode(&matcode);
  mm_set_matrix(&matcode);
  mm_set_coordinate(&matcode);
  mm_set_real(&matcode);

  //create a file and wrote into it.
    //    std::cout<<filename<<", "<<M<<", "<<N<<", "<<nz<<std::endl;
  FILE *fout;
  fout=fopen(filename,"w");
  mm_write_banner(fout, matcode);
  mm_write_mtx_crd_size(fout, M, N, nz);

  /* NOTE: matrix market files use 1-based indices, i.e. first element
      of a vector has index 1, not 0.  */

  //it is not necessary to seperate the reading and writing in 2 loops
  for (int i=0; i<nz; i++){
    fprintf(fout, "%d %d %g\n", I[i]+1, J[i]+1, val[i]);
    //originally use %10.3g, now use %g
  }
  fclose(fout);
  std::cout<<"wrote the matrix into file "<<filename<<std::endl;
  return 0;
	
}


int GF2mat_to_MM(itpp::mat G, std::string file_name){
  char temp[file_name.length()+1];
  strcpy(temp, file_name.c_str());
  return GF2mat_to_MM(G, temp);
}


int mat_to_MM(itpp::mat G, std::string file_name){
  char temp[file_name.length()+1];
  strcpy(temp, file_name.c_str());
  return mat_to_MM(G, temp);
}
