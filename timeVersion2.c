
/*
  Chloe VanCory & Kalyn Howes 
  COSC 420 
  11/ 8/2021
  Lab 4 - page rank power method

GOAL : uses MPI_File access methods to read/write the distributed “chunks” of A and x to files instead 
of using scatter/gather.
*/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include <mpi.h>

#include "matrixFunctions.c"



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
 

  // generate the numbers for matrix

 double* arrNum = malloc(rows * cols * sizeof(double));
  
  if (rank == ROOT) {
    srand(time(0));
    for (int i = 0; i < rows * cols; i++) {
      arrNum[i] = rand() % 100 + 1;
      // arrNum[i] =1;
      printf("arr[%d]=%f\n", i,arrNum[i]);
    }
  }
  MPI_Bcast(arrNum, rows * cols, MPI_DOUBLE, ROOT, world);

  int status = MPI_File_open(world,                            // comm
                             matfilename,                         // filename
                             MPI_MODE_CREATE | MPI_MODE_RDWR,  // mode
                             MPI_INFO_NULL,                    // info structure
                             &matDataFile);

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
                  rows * cols, MPI_DOUBLE, MPI_STATUS_IGNORE);
  }
  

  MPI_File_seek(matDataFile, 0, MPI_SEEK_SET);
  MPI_File_seek(matDimFile, 0, MPI_SEEK_SET);


  int totalRows, totalCols; 


  // parse the file and into a matrix of doubles
  MPI_File_read(matDimFile, 
                &totalRows, 
                1, 
                MPI_INT, 
                MPI_STATUS_IGNORE);
//   printf("Rank %d A.rows =%d\n", rank, totalRows);

  MPI_File_read(matDimFile, 
                &totalCols, 
                1, 
                MPI_INT, 
                MPI_STATUS_IGNORE);
//   printf("Rank %d A.cols =%d\n", rank, totalCols );



//   // everyone knows how many rows and cols to read in 
  SGData local_count = getSGCounts(totalRows, totalCols , worldSize);

//   // for(int i = 1; i< worldSize ;i++){
//     // printf("rank=%d offsetArr[]=%d\n",rank, local_count.cnts[rank]);
//   // }

  Matrix A;
  A.rows = local_count.cnts[rank]/totalCols;
  A.cols = totalCols ;
  A.data = malloc(A.rows * A.cols * sizeof(double));
//   printf("rank=%d rows=%D cols=%d\n", rank,A.rows,A.cols);

  int * offset = malloc(worldSize * sizeof(int));
  offset[0] = 0;
  // offset[1] = offset[0] + sizeof(double) * local_count.cnts[0];
  // offset[2] = offset[1] + sizeof(double) * local_count.cnts[1];
  for(int i = 1; i< worldSize ;i++){
    offset[i] = offset[i-1] + sizeof(double) * local_count.cnts[i-1];

  }
//   printf("rank=%d offset=%d \n", rank,offset[rank]);

    
  MPI_File_read_at(matDataFile,offset[rank], A.data,
                    local_count.cnts[rank], MPI_DOUBLE, MPI_STATUS_IGNORE);
    // if(rank ==0 )printMatrix(A);
    // if(rank ==1)printMatrix(A);

//   // rank makes 1's vector - X- and then Bcast 
  Matrix onesVector;
  onesVector.rows = totalRows;
  onesVector.cols = 1;
  onesVector.data = malloc(onesVector.rows* onesVector.cols*sizeof(double));
//   printf("rank=%d onesVectorRows=%D onesVectorCols=%d\n", rank,onesVector.rows,onesVector.cols);

  if(rank ==0){
    for(int i = 0; i < onesVector.rows * onesVector.cols ; i ++){
        onesVector.data[i] = 1;
    }

  }
  

  MPI_Bcast(onesVector.data,onesVector.rows* onesVector.cols , MPI_DOUBLE, ROOT, world);
  


  double startTime, stopTime;
  startTime = MPI_Wtime();
  double eigenvalue = powerMethod(A ,onesVector, totalRows, totalCols , THRESHOLD, e);
  if(rank ==0 )printf("EIGENVALUE= %f\n",eigenvalue );
  stopTime = MPI_Wtime();
  if(rank ==0 ) printf("MPI_Wtime measured: %2.5f seconds\n", stopTime-startTime); 
  fflush(stdout); // manually flush the buffer for safety

  MPI_File_close(&matDataFile);
  MPI_File_close(&matDimFile);

  MPI_Finalize();
  return 0;
}