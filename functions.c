/* CS 140
 * Assignment 2 : Matrix Vector Multiplication and the Power Method 
 * Group members : <Team-member-1> , <Team-member-2>
 * */

/* This file provides the placeholder function definitions, where you will be
 * implementing the assignment algorithms. You will be required to turn in 
 * only this file during the submission, where it will be compiled together
 * with our main function and tested. It is therefore required that you keep the 
 * function declaration formats unchanged.
 */

#include "powermethod.h"

// Subroutine for generating the input matrix (just one thread's part)
void generatematrix(double * mat, int size)
{
    //fill in mat with values
    //mat[i][k]
    //processor 0 gets rows 1 through n/p of A, processor 1 gets rows n/p + 1 through 2n/p, and so forth
    
    int myrank,numprocs;
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    
    for( int i = myrank*(size/numprocs)+1 -1; i < (myrank+1)(n/p); i++ ){
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
  return lambda;
}

double norm2(double *x, int size);

void matVec(double *mat, double *vec, double *local_vec, int nrows, int size);
