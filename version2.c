
/*
Chloe VanCory & Kalyn Howes
Lab 4 

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

  int totalRows, totalCols; 


  // parse the file and into a matrix of doubles
  MPI_File_read(matDimFile, 
                &totalRows, 
                1, 
                MPI_INT, 
                MPI_STATUS_IGNORE);
  // printf("Rank %d A.rows =%d\n", rank, A.rows);

  MPI_File_read(matDimFile, 
                &totalCols, 
                1, 
                MPI_INT, 
                MPI_STATUS_IGNORE);
  // printf("Rank %d A.cols =%d\n", rank, A.cols );




  // everyone knows how many rows and cols to read in 
  SGData local_count = getSGCounts(totalCols, totalCols , worldSize);

  // for(int i = 1; i< worldSize ;i++){
    // printf("rank=%d offsetArr[]=%d\n",rank, local_count.cnts[rank]);
  // }

  Matrix A;
  A.rows = local_count.cnts[rank]/totalCols;
  A.cols = totalCols ;
  A.data = malloc(A.rows * A.cols * sizeof(double));
  // printf("rank=%d rows=%D cols=%d\n", rank,A.rows,A.cols);

  int * offset = malloc(worldSize * sizeof(int));
  offset[0] = 0;
  // offset[1] = offset[0] + sizeof(double) * local_count.cnts[0];
  // offset[2] = offset[1] + sizeof(double) * local_count.cnts[1];
  for(int i = 1; i< worldSize ;i++){
    offset[i] = offset[i-1] + sizeof(double) * local_count.cnts[i-1];
  }
    
  MPI_File_read_at(matDataFile,offset[rank], A.data,
                    local_count.cnts[rank], MPI_DOUBLE, MPI_STATUS_IGNORE);


  // rank makes 1's vector - X- and then Bcast 
  Matrix onesVector;
  onesVector.rows = totalRows;
  onesVector.cols = 1;
  onesVector.data = malloc(onesVector.rows* onesVector.cols*sizeof(double));
  
  for(int i = 0; i < onesVector.rows * onesVector.cols ; i ++){
    onesVector.data[i] = 1;
  }

  MPI_Bcast(onesVector.data,onesVector.rows* onesVector.cols , MPI_DOUBLE, ROOT, world);
  
  double eigenvalue = powerMethod(A,onesVector,totalRows, totalCols , THRESHOLD,e);
  printf("EIGENVALUE= %f\n",eigenvalue );


  MPI_File_close(&matDataFile);
  MPI_File_close(&matDimFile);

  MPI_Finalize();
  return 0;
}