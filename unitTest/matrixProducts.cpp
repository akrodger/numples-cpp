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

#define MIN_ALLOWED_ARGS 5

using namespace numples;

/*
 * Function declarations go here:
 */

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
	double r1, r2, dpVal=0;//used for testing dot product.
	lapack_int m, n, p, q, i, j, k;//matrix dimensions and loop iterators.
	if(argc < MIN_ALLOWED_ARGS){
		printf("\nA small script to test the various matrix products.");
		printf("\n===================================================");
		printf("\n%s [m] [n] [p] [q]",argv[0]);
		printf("\nDoes a test of dot product, matrix multiply,\n"
				"hadamard product, and kronecker product.\n");
		exit(1);
	}
	m = (lapack_int) strtol(argv[1], &pEnd, 0);
	n = (lapack_int) strtol(argv[2], &pEnd, 0);
	p = (lapack_int) strtol(argv[3], &pEnd, 0);
	q = (lapack_int) strtol(argv[4], &pEnd, 0);
	Matrix A = Matrix(m,n);
	Matrix B = Matrix(p,q);
	Matrix C;
	Matrix K;
	Matrix u = Matrix(n,1);
	Matrix v = Matrix(n,1);
	for(j=0;j<n;j++){
		r1 = ((double)rand()) / (RAND_MAX + 1.0);
		r2 = ((double)rand()) / (RAND_MAX + 1.0);
		dpVal += r1*r2;
		u(j) = r1;
		v(j) = r2;
		for(i=0;i<m;i++){
			A(i,j) = ((double)rand()) / (RAND_MAX + 1.0);
		}
	}
	for(k=0;k<q;k++){
		for(j=0;j<p;j++){
			B(j,k) = ((double)rand()) / (RAND_MAX + 1.0);
		}
	}
	printf("Testing Dot Product\nErr = %le"
	"\nSuccess if zero.\n", ABS(u.dot(v) - dpVal));
	printf("\nTesting Matrix Multiplication.\n");
	C = A*B;
	printf("\nIf no error + crash, successful.\n");
	printf("\nTesting Hadamard Product.\n");
	C = A.hada(B);
	printf("\nIf no error + crash, successful.\n");
	printf("\nTesting Kronecker Product.\n");
	K = A.kron(B);
	printf("\nIf no error + crash, successful.\n");
	printf("\n%s\n\n", argv[0]);
	return 0;
}
 
