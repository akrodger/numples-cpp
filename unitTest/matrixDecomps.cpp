/*
 * C Driver file Template by Bram Rodgers.
 * Original Draft Dated: 25, Feb 2018
 */

/*
 * Macros and Includes go here: (Some common ones listed)
 */
#include<stdio.h>
#include<stdlib.h>
#include"../src/numples.h"
#ifndef NULL
	#define NULL 0
#endif
#define MIN_ALLOWED_ARGS 4
/*
 * Function declarations go here:
 */
using namespace numples;
/*
 * Template main:
 */
int main(int argc, char** argv){
	// Variable Declarations:
	// int, double, etc
	// Pre Built argv parsing: (argv[0] is executable title. e.g. "./a.out")
	char* pEnd = NULL; //Points to end of a parsed string
	/*
	 * (intger) = atoi(argv[index]);
	 * (double) = strtod(argv[index], &pEnd);
	 * (long int) = strtol(argv[index], &pEnd, base_val)
	 */
	srand(time(NULL));
	lapack_int m, n, numSV;//matrix dimensions and loop iterators.
	if(argc < MIN_ALLOWED_ARGS){
		printf("\nA small script to test the various matrix products.");
		printf("\n===================================================");
		printf("\n%s [m] [n] [numSV]",argv[0]);
		printf("\nDoes a test of SVD, obtaining only numSV singular vectors,\n"
				"and Schur Decomposition.\n");
		exit(1);
	}
	m = (lapack_int) strtol(argv[1], &pEnd, 0);
	n = (lapack_int) strtol(argv[2], &pEnd, 0);
	numSV = (lapack_int) strtol(argv[3], &pEnd, 0);
	Matrix A = Matrix(m , n, 'U');//generate uniform random matrix.
	Matrix R, Z;//Schur decomposition matrices
	Matrix U, sig, Vt;//SVD matrices
	Matrix Resid;
	double residNorm, traceVal;
	A.svd(U, sig, Vt);
	Resid = A - U * sig.diag() * Vt;
	residNorm = Resid.norm2();
	printf("\nSVD Residual: %le\n", residNorm);
	Z = A.getLeftSV(numSV);
	traceVal = (Z.transp() * U).trace();
	residNorm = numSV - ABS(traceVal);
	printf("\n%le\n",residNorm);
	if(residNorm <= (double)numSV){
		printf("\nGet Left Singular Vectors passed within eps_mach: 1\n");
	}else{
		printf("\nGet Left Singular Vectors failed within eps_mach: 0\n");
	}
	A.schur(Z,R);
	Resid = A*Z - Z*R;
	residNorm = Resid.norm2();
	printf("\nSchur Residual: %le\n", residNorm);
	printf("\n%s\n\n", argv[0]);
	return 0;
}
 
