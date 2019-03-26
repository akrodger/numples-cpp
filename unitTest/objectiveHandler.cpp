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
#include"../src/numples-numOpt.h"
#ifndef NULL
	#define NULL 0
#endif

#define TEST_DIM 5

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
	Matrix u = Matrix(TEST_DIM,1,'U');
	Matrix v = Matrix(TEST_DIM,1,'U');
	Matrix H, Hv, g, r;
	u = u*5; //make u a little bigger
	Objective J = Objective(u);
	g = J.grad();
	r = g - u;
	printf("\nApprox val of grad resid    : %le\n", sqrt(r.dot(r)));
	printf("\nStored Value of Objective   : %le\n", J.getValue());
	printf("\nComputed Value of Objective : %le\n", J(u));
	r = u + g*(-1e-3);
	J.update(r);
	printf("\nSmaller Value of Objective  : %le\n", J.getValue());
	printf("\nIs it actually smaller?     : %d\n", J.getValue() < J(u));
	printf("\nApproximate Full Hessian    :\n");
	H = J.hess();
	H.print();
	return 0;
}
