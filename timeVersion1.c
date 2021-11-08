
/*

  Chloe VanCory & Kalyn Howes 
  COSC 420 
  11/ 8/2021
  Lab 4 - page rank power method

GOAL : the root processor reads the entire A and x data from a file and scatters them appropriately to 
perform the matrix multiplication with your existing library. Then use other collectives as appropriate to 
perform the update and normalization in parallel.

OUTPUT

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
  int rows = atoi(argv[1]);
  int cols = atoi(argv[2]);
  int THRESHOLD = atoi(argv[3]);
  double e = 10E-16;


  if(argc != 4){
    //   MPI_Barrier(world);
      printf("NOT ENOUGH ARGS PASSED FOR DIMENSIONS");
      MPI_Finalize();
      return 1;
  }

  
  // generate the numbers
  double* arrNum = malloc(rows * cols * sizeof(double));
  
  if (rank == ROOT) {
    srand(time(0));
    for (int i = 0; i < rows * cols; i++) {
      arrNum[i] = rand() % 100 + 1;
      // arrNum[i] =1;
    //   printf("arrNum[i]=%f\n", arrNum[i]);
    }
  }
  MPI_Bcast(arrNum, rows * cols, MPI_DOUBLE, ROOT, world);
  

  // // open the file
  int status = MPI_File_open(world,                            // comm
                             matfilename,                         // filename
                             MPI_MODE_CREATE | MPI_MODE_RDWR,  // mode
                             MPI_INFO_NULL,                    // info structure
                             &matDataFile);

  if (status == -1) {
    MPI_Finalize();
    puts("ERROR opening file");
     MPI_Finalize();
    return 1;
  }

  status = MPI_File_open(world,                            // comm
      dimfilename,                         // filename
      MPI_MODE_CREATE | MPI_MODE_RDWR,  // mode
      MPI_INFO_NULL,                    // info structure
      &matDimFile);

  
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

  
  
  Matrix A;

  // parse the file and into a matrix of doubles

  MPI_File_read(matDimFile, 
                &A.rows, 
                1, 
                MPI_INT, 
                MPI_STATUS_IGNORE);
  printf("Rank %d A.rows =%d\n", rank, A.rows);

  MPI_File_read(matDimFile, 
                &A.cols, 
                1, 
                MPI_INT, 
                MPI_STATUS_IGNORE);
  printf("Rank %d A.cols =%d\n", rank, A.cols );

  A.data = malloc(A.rows * A.cols * sizeof(double));

  MPI_File_read(matDataFile, 
                A.data, 
                A.rows*A.cols, 
                MPI_DOUBLE, 
                MPI_STATUS_IGNORE);


  if(rank ==ROOT){
    printMatrix(A);
  }


  SGData scatter = getSGCounts(A.rows, A.cols , worldSize);

  Matrix localMatrix;
  localMatrix.rows =  scatter.cnts[rank] / A.rows;
  localMatrix.cols = A.cols;
  localMatrix.data = malloc(localMatrix.rows  * localMatrix.cols * sizeof(double));
//   printf("RANK =%d localMatrix rows= %d cols=%d\n", rank,localMatrix.rows,localMatrix.cols);

  MPI_Scatterv(
            A.data,                // sendbuf
            scatter.cnts,        // sendcnts
            scatter.displs,      // displacements
            MPI_DOUBLE,               // datatype
            localMatrix.data,              // recvbuf
            scatter.cnts[rank],  // recvcnt
            MPI_DOUBLE, ROOT, world);


  // char *arrbuf = bufArr(localMatrix.data, scatter.cnts[rank]);
  // printf("Rank %d received %s\n", rank, arrbuf);  

  // everyone declares a  X vector 
  Matrix X; 
  X.rows = localMatrix.cols;
  X.cols = 1;
  X.data= malloc(X.rows * X.cols * sizeof(double));
//   printf("RANK =%d localX rows= %d cols=%d\n", rank,X.rows,X.cols);


  /* set X to the 1's vector */
  for(int i =0 ;i < X.rows* X.cols; i++) X.data[i] = 1;
//   char *arrbuf2 = bufArr(X.data, X.rows* X.cols);
//   printf("Rank %d received %s\n", rank, arrbuf2);  


  double startTime, stopTime;
  startTime = MPI_Wtime();
  double eigenvalue = powerMethod(localMatrix,X,A.rows, A.cols, THRESHOLD, e);
  stopTime = MPI_Wtime();
  if(rank ==0 ) printf("MPI_Wtime measured: %2.5f seconds\n", stopTime-startTime); 
  fflush(stdout); // manually flush the buffer for safety

  if(rank ==0 )printf("EIGENVALUE= %f\n",eigenvalue );
  


  MPI_File_close(&matDataFile);
  MPI_File_close(&matDimFile);
  
  MPI_Finalize();
  return 0;
}

