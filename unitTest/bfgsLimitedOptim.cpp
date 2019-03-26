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

using namespace numples;

/*
 * Function declarations go here:
 */

class QuarticObjective : public Objective{
private:
	double hessCurve;
public:
	//Default constructor sets epsi to machine accuracy, val=0, x to zero.
	QuarticObjective(double hC) : Objective() {
		this->hessCurve = hC;
	};
	//A constructor which allows you to set initial guess yourself.
	//val is set to this->operator()(x).
	QuarticObjective(double hC, const Matrix& x0){
		this->hessCurve = hC;
		this->Objective::update(x0);
	}
	//A constructor which allows you to set epsi and initial guess yourself.
	//val is set to this->operator()(x). If e is negative, it defaults to
	//machine accuracy.
	QuarticObjective(double hC, double e, const Matrix& x0) : Objective(e) {
		this->hessCurve = hC;
		this->Objective::update(x0);
	}
	virtual double operator()(const Matrix& y){
		double c = this->hessCurve;
		double term1 = (c*y.get(0) - 2);
		double term2 = term1*term1;
		double term3 = (y.get(1) + 1);
		term1 = term2*term2;
		term2 = term2*y.get(1)*y.get(1);
		term3 = term3*term3;
		return term1 + term2 + term3;
	}
};


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
	srand(time(0));
	if(argc < 2){
		printf("\nA small script to test Wolfe Line Search.");
		printf("\nThis is accomplished by solving a 2D gradient descent.");
		printf("\n=======================================================");
		printf("\n%s [hessCurve]",argv[0]);
		printf("\nhessCurve is type double. Bigger means harder to solve.\n\n");
		exit(1);
	}
	int stopSignal = 0;
	lapack_int numIter = 0;
	double hC = strtod(argv[1], &pEnd);
	Matrix initGuess = Matrix(2,1,'U');
	printf("\nInitial Guess:\n");
	initGuess.print();
	Matrix xOpt;
	Matrix currGrad;
	QuarticObjective J = QuarticObjective(hC, 1e-10, initGuess);
	printf("\nInitial Objective value: %le\n", J.getValue());
	//currGrad = J.grad();
	//currGrad.print();
	numIter = bfgsLimitedMem(	J,
								1e-4,
								0.9,
								1e-4,
								stopSignal,
								10);
	printf("\nbfgs Objective Value: %le\n", J.getValue());
	J.getArg().print();
	J.grad().print();
	printf(	"\nStop signal was : %d\n"
			  "Iterations was  : %ld\n", stopSignal, numIter);
	printf("\n%s\n\n", argv[0]);
	return 0;
}
 
