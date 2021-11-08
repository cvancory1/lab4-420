
/*

Chloe VanCory & Kalyn Howes
Lab 4 - 

GOAL : the root processor reads the entire A and x data from a file and scatters them appropriately to 
perform the matrix multiplication with your existing library. Then use other collectives as appropriate to 
perform the update and normalization in parallel.

*/

#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>


#include "matrixFunctions.c"


void printbuffer(double * buf, int size){
    for(int i =0; i < size ;i++){
        printf("arr=%f\n", buf[i]);
    }

}

// void printMatrxD(Matrix A) {
//   for (int i = 0; i < A.rows; i++) {
//     for (int j = 0; j < A.cols; j++) {
//       printf("%2f ", ACCESS(A, i, j));
//     }
//     printf("\n");
//   }
// }

Matrix copyMatrix(Matrix copyMatrix){
    Matrix newMatrix;
    newMatrix.rows = copyMatrix.rows ;
    newMatrix.cols = copyMatrix.cols;
    newMatrix.data = malloc(newMatrix.rows *  newMatrix.cols  * sizeof(double));

    for(int i =0; i < copyMatrix.rows *  copyMatrix.cols ;i++ ){
      newMatrix.data[i] = copyMatrix.data[i];
    }

    return newMatrix;

}

int MIN(int a , int b){
  if(a<b){
    return a;
  }else{
    return b;
  }
}



int main(int argc, char** argv) {
  
  MPI_Init(&argc, &argv);

  world = MPI_COMM_WORLD;

  // passing the container for the result in second param
  MPI_Comm_size(world, &worldSize);
  MPI_Comm_rank(world, &rank);
  MPI_Get_processor_name(name, &nameLen);

  MPI_File matDataFile, matDimFile;

  char* matfilename = "matData.txt";
  char* dimfilename = "matDims.txt";
  int rows = 5;
  int cols = 5;
  int THRESHOLD = 20;
  double e = 10E-16;
  double eigenvalue ;


  /*
  // generate the numbers
  double* arrNum = malloc(rows * cols * sizeof(double));
  
  if (rank == ROOT) {
    srand(time(0));
    for (int i = 0; i < rows * cols; i++) {
      arrNum[i] = rand() % 100 + 1;
      // arrNum[i] =1;
      printf("arrNum[i]=%f\n", arrNum[i]);
    }
  }
  MPI_Bcast(arrNum, rows * cols, MPI_DOUBLE, ROOT, world);
  */

  // // open the file
  int status = MPI_File_open(world,                            // comm
                             matfilename,                         // filename
                             MPI_MODE_CREATE | MPI_MODE_RDWR,  // mode
                             MPI_INFO_NULL,                    // info structure
                             &matDataFile);

  if (status == -1) {
    MPI_Finalize();
    puts("ERROR opening file");
    return 1;
  }

  status = MPI_File_open(world,                            // comm
      dimfilename,                         // filename
      MPI_MODE_CREATE | MPI_MODE_RDWR,  // mode
      MPI_INFO_NULL,                    // info structure
      &matDimFile);

  /*
  // // write to the newfile
  MPI_Barrier(world);
  if( rank == 0 ){
  int offset = 0;
  MPI_File_write(matDimFile,
                    // 1,  // offset
                    &rows,   // buf
                    1, MPI_INT, MPI_STATUS_IGNORE);

  offset += sizeof(int);
   MPI_File_write(matDimFile,
                    // 1,  // offset
                    &cols,   // buf
                    1, MPI_INT, MPI_STATUS_IGNORE);

  offset += sizeof(int);
    MPI_File_write(matDataFile,
                  // 1,  // offset
                  arrNum,   // buf
                  rows *cols, MPI_DOUBLE, MPI_STATUS_IGNORE);
  }
  

  MPI_File_seek(matDataFile, 0, MPI_SEEK_SET);
  MPI_File_seek(matDimFile, 0, MPI_SEEK_SET);

  */
  
  // Matrix A;


  // // parse the file and into a matrix of doubles

  // MPI_File_read(matDimFile, 
  //               &A.rows, 
  //               1, 
  //               MPI_INT, 
  //               MPI_STATUS_IGNORE);
  // // printf("Rank %d A.rows =%d\n", rank, A.rows);

  // MPI_File_read(matDimFile, 
  //               &A.cols, 
  //               1, 
  //               MPI_INT, 
  //               MPI_STATUS_IGNORE);
  // // printf("Rank %d A.cols =%d\n", rank, A.cols );

  // A.data = malloc(A.rows * A.cols * sizeof(double));

  // MPI_File_read(matDataFile, 
  //               A.data, 
  //               A.rows*A.cols, 
  //               MPI_DOUBLE, 
  //               MPI_STATUS_IGNORE);


  // if(rank ==ROOT){
  //   printMatrix(A);
  // }


  // SGData scatter = getSGCounts(A.rows, A.cols , worldSize);

  // Matrix localMatrix;
  // localMatrix.rows =  scatter.cnts[rank] / A.rows;
  // localMatrix.cols = A.cols;
  // localMatrix.data = malloc(localMatrix.rows  * localMatrix.cols * sizeof(double));
  // // printf("RANK =%d localMatrix rows= %d cols=%d\n", rank,localMatrix.rows,localMatrix.cols);

  // MPI_Scatterv(
  //           A.data,                // sendbuf
  //           scatter.cnts,        // sendcnts
  //           scatter.displs,      // displacements
  //           MPI_DOUBLE,               // datatype
  //           localMatrix.data,              // recvbuf
  //           scatter.cnts[rank],  // recvcnt
  //           MPI_DOUBLE, ROOT, world);


  // char *arrbuf = bufArr(localMatrix.data, scatter.cnts[rank]);
  // printf("Rank %d received %s\n", rank, arrbuf);  

  // everyone declares a  X vector 
  // Matrix X; 
  // X.rows = localMatrix.cols;
  // X.cols = 1;
  // X.data= malloc(X.rows * X.cols * sizeof(double));
  // // printf("RANK =%d localX rows= %d cols=%d\n", rank,X.rows,X.cols);


  // /* set X to the 1's vector */
  // for(int i =0 ;i < X.rows* X.cols; i++) X.data[i] = 1;
  // char *arrbuf2 = bufArr(X.data, X.rows* X.cols);
  // // printf("Rank %d received %s\n", rank, arrbuf2);  

 
  //  eigenvalue = powerMethod(localMatrix,X, A.rows, A.cols,THRESHOLD,e);
  // if(rank ==0 )printf("EIGENVALUE= %f\n",eigenvalue );







  /* BELOW IS HARD CODED EXAMPLES FOR TESTING */ 
  /* testing kalyns matrix */ 
  double arr2[] = {0,1,-2,-3};
  Matrix Kalyn;
  Kalyn.rows =2;
  Kalyn.cols =2;
  Kalyn.data = malloc(Kalyn.rows* Kalyn.cols*sizeof(double));
  
  for(int i = 0; i < Kalyn.rows * Kalyn.cols ; i ++){
    Kalyn.data[i] = arr2[i];
  }


  SGData K_counts = getSGCounts(Kalyn.rows, Kalyn.cols , worldSize);

  // printf("counsts= %d displs=%d\n",K_counts.cnts[rank], K_counts.displs[rank] );
  Matrix onesVector;
  onesVector.rows = Kalyn.rows;
  onesVector.cols = 1;
  onesVector.data = malloc(onesVector.rows* onesVector.cols*sizeof(double));
  
  for(int i = 0; i < onesVector.rows * onesVector.cols ; i ++){
    onesVector.data[i] = 1;
  }

  MPI_Bcast(onesVector.data,onesVector.rows* onesVector.cols , MPI_DOUBLE, ROOT, world);
  // matrixMultiplication(Kalyn,onesVector);
  
  Matrix localK;
  localK.rows = K_counts.cnts[rank] / Kalyn.rows;
  localK.cols = Kalyn.cols;
  localK.data = malloc(Kalyn.rows* Kalyn.cols* sizeof(double));

  MPI_Scatterv(
          Kalyn.data,                // sendbuf
          K_counts.cnts,        // sendcnts
          K_counts.displs,      // displacements
          MPI_DOUBLE,               // datatype
          localK.data,              // recvbuf
          K_counts.cnts[rank],  // recvcnt
          MPI_DOUBLE, ROOT, world);
  // char *arrbuf2 = bufArr(localK.data, localK.rows* localK.cols);
  // printf("Rank %d received %s\n", rank, arrbuf2);  
 
  eigenvalue =  powerMethod(localK,onesVector, Kalyn.rows, Kalyn.cols,THRESHOLD,e);
  if(rank ==0 )printf("EIGENVALUE= %f\n",eigenvalue );




  double arr3[] = {
                  6,3,6,5,2,
                  9,8,8,6,6,
                  8,3,3,8,3,
                  4,1,3,8,9,
                  10,8,4,5,1
  };
  Matrix Bevan;
  Bevan.rows =5;
  Bevan.cols =5;
  Bevan.data = malloc(Bevan.rows* Bevan.cols* sizeof(double));
   for(int i = 0; i < Bevan.rows * Bevan.cols ; i ++){
    Bevan.data[i] = arr3[i];
  }

  SGData B_counts = getSGCounts(Bevan.rows, Bevan.cols , worldSize);

  Matrix localB;
  localB.rows = B_counts.cnts[rank] / Bevan.rows;
  localB.cols = Bevan.cols;
  localB.data = malloc(Bevan.rows* Bevan.cols* sizeof(double));


  Matrix onesVector2;
  onesVector2.rows = Bevan.rows;
  onesVector2.cols = 1;
  onesVector2.data = malloc(onesVector2.rows* onesVector2.cols*sizeof(double));
  
  for(int i = 0; i < onesVector2.rows * onesVector2.cols ; i ++){
    onesVector2.data[i] = 1;
  }


 MPI_Scatterv(
          Bevan.data,                // sendbuf
          B_counts.cnts,        // sendcnts
          B_counts.displs,      // displacements
          MPI_DOUBLE,               // datatype
          localB.data,              // recvbuf
          B_counts.cnts[rank],  // recvcnt
          MPI_DOUBLE, ROOT, world);
 
  eigenvalue = powerMethod(localB,onesVector2, Bevan.rows, Bevan.cols, THRESHOLD,e);
  if(rank ==0 )printf("EIGENVALUE= %f\n",eigenvalue );



/* KALYN 3x3 trace */ 
  double arr4[] = {
                 2,8,10,
                 8,3,4,
                 10,4,7
  };
  Matrix test3;
  test3.rows =3;
  test3.cols =3;
  test3.data = malloc(test3.rows* test3.cols* sizeof(double));
   for(int i = 0; i < test3.rows * test3.cols ; i ++){
    test3.data[i] = arr4[i];
  }

  if(rank==0) printMatrix(test3);

  SGData test3_counts = getSGCounts(test3.rows, test3.cols , worldSize);

  Matrix localtest3;
  localtest3.rows = test3_counts.cnts[rank] / test3.rows;
  localtest3.cols = test3.cols;
  localtest3.data = malloc(test3.rows* test3.cols* sizeof(double));
  // printf("RANK =%d localtest3 rows= %d cols=%d\n", rank,localtest3.rows,localtest3.cols);


  Matrix onesVector3;
  onesVector3.rows = test3.rows;
  onesVector3.cols = 1;
  onesVector3.data = malloc(onesVector3.rows* onesVector3.cols*sizeof(double));
  
  for(int i = 0; i < onesVector3.rows * onesVector3.cols ; i ++){
    onesVector3.data[i] = 1;
  }


 MPI_Scatterv(
          test3.data,                // sendbuf
          test3_counts.cnts,        // sendcnts
          test3_counts.displs,      // displacements
          MPI_DOUBLE,               // datatype
          localtest3.data,              // recvbuf
          test3_counts.cnts[rank],  // recvcnt
          MPI_DOUBLE, ROOT, world);
 
  // if(rank ==0) printMatrix(localtest3);
  // puts("");
  // if(rank ==1) printMatrix(localtest3);

  eigenvalue = powerMethod(localtest3,onesVector3, test3.rows, test3.cols, 3,THRESHOLD);
  printf("EIGENVALUE= %f\n",eigenvalue );




  


  MPI_File_close(&matDataFile);
  MPI_File_close(&matDimFile);
  
  MPI_Finalize();
  return 0;
}










