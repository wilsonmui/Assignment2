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

// Subroutine for generating the input matrix (just one thread's part)
void generatematrix(double * mat, int size)
{
    //fill in mat with values
    //mat[i][k]
    //processor 0 gets rows 1 through n/p of A, processor 1 gets rows n/p + 1 through 2n/p, and so forth
    
    int myrank,numprocs;
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    
    for( int i = myrank*(size/numprocs); i < (myrank+1)(n/p); i++ ){
        for( int k = 0; k < i; k++){
            mat[i][k] = i;
        }
    }
}

// Subroutine to generate a vector
void generatevec(double * x,int size)
{
    int myrank;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    
    for (int i = 0; i < size; i++){
        x[i] = 1;
    }
    
}

// Subroutine for the power method, to return the spectral radius
double powerMethod(double * mat, double * x, int size, int iter)
{
    int myrank, numprocs, product_vector[size], result_vec, lambda;
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    
    //each process should broadcast and gather
    /*
     int arr[size];
     //mult
     for(int i = 0; i < numprocs; i++){
        MPI_Bcast(arr[i * row_per_proc], row_per_proc, MPI_INT, i, MPI_COMM_WORLD);
     }
     //norm
     */
    
    int num_rows = size/numprocs;
    
    //initialize result_vec to all 1's
    generatevec(result_vec, size);
    for(int k = 0; k < iter; k++){
        lambda = norm2(result_vec, size);
        
        for (int i = 0; i < size; i++){
            x[i] = x[i]/lambda;
        }
        
        //make local vec = matrix * local vec
        matVec(mat, x, x, num_rows, size);
        
        //make all processes' local vectors contain the full product vector
        /*
        for(int i = 0; i < numprocs; i++){
            for(int j = myrank*(size/numprocs); j < (myrank+1)(n/p); j++){
                MPI_Bcast(&x[j], 1, MPI_INT, i, MPI_COMM_WORLD);
            }
        }
         */
        for(int i = 0; i < numprocs; i++){
            MPI_Bcast(x[i * row_per_proc], row_per_proc, MPI_INT, i, MPI_COMM_WORLD);
        }
    }
    
}

//compute the 2-norm (length) of a given vector
double norm2(double *x, int size){
    int sum = 0;
    for(int i = 0; i < size; i++){
        r += x[i]*x[i];
    }
    return sqrt(r);
}

//multiply matrix by a vector.
//should result in vec being result
void matVec(double *mat, double *vec, double *local_vec, int nrows, int size){
    int myrank,numprocs;
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    
    for(int j = myrank*(size/numprocs); j < (myrank+1)(n/p); j++){
        for (int k = 0; k < nrows; k++){
            for (int i = 0; i < size; i++){
                vec[j] += mat[k][i] * local_vec[i];
            }
        }
    }
}


















