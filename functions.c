/* CS 140
 * Assignment 2 : Matrix Vector Multiplication and the Power Method 
 * Group members : Wilson Mui , Karl Wang
 * */

/* This file provides the placeholder function definitions, where you will be
 * implementing the assignment algorithms. You will be required to turn in 
 * only this file during the submission, where it will be compiled together
 * with our main function and tested. It is therefore required that you keep the 
 * function declaration formats unchanged.
 */

#include "powermethod.h"
#include <math.h>

const int DEBUG = 1;

// Subroutine for generating the input matrix (just one thread's part)
void generatematrix(double * mat, int size)
{
    //fill in mat with values
    //mat[i][k]
    //processor 0 gets rows 1 through n/p of A, processor 1 gets rows n/p + 1 through 2n/p, and so forth
    
    int myrank, numprocs, row, col;
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    int num_rows = size/numprocs;
    
    if(DEBUG) { printf("matrix of process %d/%d:\n", myrank, numprocs); }
    for( row = 0; row < num_rows; row++ ){
        int realrow = myrank*num_rows + row;
        for( col = 0; col < size; col++ ){
            //mat[row*size + col] = ( (row+11)*(col+13)*(myrank*17) ) % 1231; //persudo random by prime
            mat[row*size + col] = col <= realrow ? realrow + 1 : 0; //persudo random by prime
            if(DEBUG) { printf("%f ", mat[row*size + col]); }
        }
        if(DEBUG) { printf("\n"); }
    }
}

// Subroutine to generate a vector
void generatevec(double * x, int size)
{
    int myrank, i;
    
    for( i = 0; i < size; i++){
        x[i] = 1;
    }
    
}

// Subroutine for the power method, to return the spectral radius
double powerMethod(double * mat, double * x, int size, int iter)
{
    int myrank, numprocs;
    double lambda;
    
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    int num_rows = size/numprocs;
    double * local_vec = calloc(num_rows, sizeof(double));
    
    //each process should broadcast and gather
    /*
     int arr[size];
     //mult
     for( i = 0; i < numprocs; i++){
        MPI_Bcast(arr[i * row_per_proc], row_per_proc, MPI_INT, i, MPI_COMM_WORLD);
     }
     //norm
     */
    
    
    int i, k;
    for( k = 0; k < iter; k++){
        lambda = norm2(x, size);
        
        for( i = 0; i < size; i++){
            x[i] = x[i]/lambda;
        }

        if( DEBUG && myrank == 0 ) {
            printf("[");
            for( i = 0; i < size; i++){
                if(DEBUG) { printf("%f ", x[i]); }
            }
            printf("] lambda: %f\n", lambda);
        }
        
        //clear local_vec
        memset(local_vec, 0, sizeof(double) * num_rows);
        //make local_vec = matrix * x
        matVec(mat, x, local_vec, num_rows, size);
        //copy values to x
        memcpy(&x[num_rows*myrank], local_vec, sizeof(double) * num_rows);
        
        //make all processes' local vectors contain the full product vector
        /*
        for( i = 0; i < numprocs; i++){
            for( j = myrank*(size/numprocs); j < (myrank+1)(n/p); j++){
                MPI_Bcast(&x[j], 1, MPI_INT, i, MPI_COMM_WORLD);
            }
        }
         */
        for( i = 0; i < numprocs; i++){
            MPI_Bcast(&x[i * num_rows], num_rows, MPI_INT, i, MPI_COMM_WORLD);
        }
    }
    
    return lambda;
    
}

//compute the 2-norm (length) of a given vector
double norm2(double *x, int size){
    int sum = 0, i;
    for( i = 0; i < size; i++){
        sum += x[i]*x[i];
    }
    return sqrt(sum);
}

//multiply matrix by a vector.
//precondition: local_vec is clean
//input: mat, vec
//output: local_vec 
void matVec(double *mat, double *vec, double *local_vec, int nrows, int size){
    int row, col;
    for( row = 0; row < nrows; row++){
        for( col = 0; col < size; col++){
            local_vec[row] += mat[row*size + col] * vec[col];
        }
    }
}
