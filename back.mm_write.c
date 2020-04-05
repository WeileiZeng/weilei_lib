/* Weilei Zeng
*   revised from Matrix Market I/O example program
*   (See http://math.nist.gov/MatrixMarket for details.)
*/
//recomend file format for Matrix Market file is .mm
//this file try to solve the segmentation fault by not defining I, J, val. which has been solved by using the pointer in C. May 2018
//trouble shooting: segmentation fault:  directory not exist; wrong folder name
#include <stdio.h>
#include <stdlib.h>
#include "mmio.h"
#include <fstream> //ofstream ifstream
#include <itpp/itbase.h>
using namespace itpp;

using namespace std;

//int GF2mat_to_MM(GF2mat G, char* file_name="mm_temp.dat")
int GF2mat_to_MM(GF2mat G, char* file_name)
{
  int M=G.rows(),N=G.cols(),nt=M*N;//nt is the total number of elements, will find nz later. nz is the number of non-zero elements
  //make nt smaller; the size of int[] should be less than about 2000000+, which is set by the hardware indtead of C++; otherwise it creat segmentation fault
  //this could be fixed by seperate int[] into smaller ones, but in this case the matrix is already too bigfor any calculation
  // cout<<"density of GF2mat: "<<G.density()<<endl;
  nt=floor(1+nt*(G.density()+0.001)  );
  cout<<"debug 0"<<endl;
  MM_typecode matcode;
  cout<<"debug 0.5"<<endl;
  //  int I[nt],J[nt];//the location of nonzero elements
  cout<<"debug 0.6"<<endl;
  //usually nt is much larger than nz, where nz/nt is the density of the code. There is some wasted memory here.
  double val[nt];//have to use double to save the binary value in order to match the format of the file, didn't affect the result so far
  cout<<"debug 1"<<endl;

    //create a file and wrote into it.
  FILE *ftemp;//a temp file to save elements first, then add the title
  char file_name_temp[255];
  sprintf(file_name_temp,"%s_temp",file_name);
  ftemp=fopen(file_name_temp,"w");

  int pos=0;//current position, which accumelate to nz
  for(int m=0;m<M;m++){
    //      cout<<"\t m="<<m<<",";
      for(int n=0;n<N;n++){
	  if(G.get(m,n)){//if it is nonzero
	    //	    I[pos]=m;
	    //J[pos]=n;
	    //val[pos]=1;//we know it is one. so dont use G.get(m,n)
	    pos++;
	    //    for (int i=0; i<nz; i++)
	    fprintf(ftemp, "%d %d %10.3g\n", m+1, n+1, 1.0);
	  }
      }
  }
  fclose(ftemp);
  ifstream infile( file_name_temp);
  cout<<"debug 2"<<endl;
  int nz=pos;//number of nonzero elements
    //cout<<"number of ones in the matrix is nz= "<<nz<<endl;
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

  string line;
  //  char * c=line;
  while (getline(infile,line)){
    fprintf(fout,"%s\n",line.c_str());
  }
  infile.close();
    //it is not necessary to seperate the reading and writing in 2 loops
      //    for (int i=0; i<nz; i++)
      //        fprintf(fout, "%d %d %10.3g\n", I[i]+1, J[i]+1, val[i]);

    fclose(fout);
    cout<<"wrote the matrix (density:"<<G.density()<<") into file "<<file_name<<endl;
    remove(file_name_temp);

    //test
    /*
    ifstream printfile(file_name);
    while (getline(printfile,line)){
      cout<<line<<endl;
      } */   
    return 0;
	
}

//int mat_to_MM(mat G, char* filename="mm_temp.dat")
int mat_to_MM(mat G, char * filename){//use this when G is not a binary matrix
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
    //    cout<<filename<<", "<<M<<", "<<N<<", "<<nz<<endl;
  FILE *fout;
  fout=fopen(filename,"w");
  mm_write_banner(fout, matcode);
  mm_write_mtx_crd_size(fout, M, N, nz);

  /* NOTE: matrix market files use 1-based indices, i.e. first element
      of a vector has index 1, not 0.  */

  //it is not necessary to seperate the reading and writing in 2 loops
  for (int i=0; i<nz; i++)
    fprintf(fout, "%d %d %g\n", I[i]+1, J[i]+1, val[i]);
  //originally use %10.3g, now use %g
  fclose(fout);
  cout<<"wrote the matrix into file "<<filename<<endl;
  return 0;
	
}

