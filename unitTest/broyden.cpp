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

class NonlinSys : public Residual{
private:
public:
	//Default constructor sets epsi to machine accuracy, val=0, x to zero.
	NonlinSys() : Residual() {};
	//A constructor which allows you to set initial guess yourself.
	//val is set to this->operator()(x).
	NonlinSys(const Matrix& x0) : Residual(x0){
	};
	//A constructor which allows you to set epsi and initial guess yourself.
	//val is set to this->operator()(x). If e is negative, it defaults to
	//machine accuracy.
	NonlinSys(double e, const Matrix& x0) : Residual(e) {
		this->Objective::update(x0);
	};
	Matrix resid(const Matrix& y){
    Matrix r = y;
    r(0) = 2.0*(y.get(0)*y.get(0) - y.get(1))*2.0*y.get(0) + 2.0*(y.get(0)-1.0);
    r(1) = 2.0*(y.get(1) - y.get(0)*y.get(0));
    return r;
  };
};



/*
 * Template main:
 */
int main(int argc, char** argv){
  // Variable Declarations:
  // int, double, etc
  // Pre Built argv parsing: (argv[0] is executable title. e.g. "./a.out")
  /*
   * (intger) = atoi(argv[index]);
   * (double) = strtod(argv[index], &pEnd);
   * (long int) = strtol(argv[index], &pEnd, base_val)
   */
  srand(time(0));
  //if(argc < 2){
	  printf("\nA small script to test Broyden Line Search.");
	  printf("\nThis is accomplished by solving a 2D nonlinear system.");
	  printf("\n=======================================================");
	  printf("\n%s\n\n",argv[0]);
	  //exit(1);
 // }
  int stopSignal = 0;
  lapack_int numIter = 0;
  Matrix initGuess = Matrix(2,1,'U');
  printf("\nInitial Guess:\n");
  initGuess.print();
  Matrix xOpt;
  Matrix currResid;
  NonlinSys R = NonlinSys(1e-10, initGuess);
  printf("\nInitial Merit value: %le\n", R.getValue());
  //currGrad = J.grad();
  //currGrad.print();
  numIter = broydenLineSearch(R,
                              1e-4,
                              0.9,
                              1e-4,
                              stopSignal);//,
                              //10);
	printf("\nBroyden Merit Value: %le\n", R.getValue());
	R.x.print();
	R.residMatrix.print();
	printf(	"\nStop signal was : %d\n"
			  "Iterations was  : %ld\n", stopSignal, numIter);
	printf("\n%s\n\n", argv[0]);
	return 0;
}
 
