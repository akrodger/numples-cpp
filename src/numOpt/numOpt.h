#ifndef NUMPLES_OPTI_H
#define NUMPLES_OPTI_H

/*
 * C Header file Template by Bram Rodgers.
 * Original Draft Dated: 25, Feb 2018
 */

/*
 * Header File Body:
 */

/*
 * Macros and Includes go here.
 */
#include"../numples.h"
#ifndef NULL
	#define NULL 0
#endif

namespace numples{
/*
 * Object and Struct Definitions:
 */

//A generic Objective class. This class can be overloaded and given
//data and a custom implementation of the evaluation operator ().
class Objective{
private:
public:
	//A small positive number. If the norm of the Objective's gradient
	//is small than this number, the gradient is treated as being zero.
	double epsi;
	//The current point in a minimization algorithm is stored.
	Matrix x;
	//The current value of the objective is stored, to ensure we don't
	//evaluate the objective more than we need to.
	double val;
	//an int for storing if the gradient has already been computed or not.
	//this is zero by default and after an update. 
	int gradAvailable;
	//A matrix for storing the gradient in the cases where is has already
	//been computed. (such as with a call to grad)
	Matrix gradMatrix;
	//Default constructor sets epsi to machine accuracy, val=0, x to zero.
	Objective();
	//A constructor which allows you to set epsi
	Objective(double e);
	//A constructor which allows you to set initial guess yourself.
	//val is set to this->operator()(x).
	Objective(const Matrix& x0);
	//A constructor which allows you to set epsi and initial guess yourself.
	//val is set to this->operator()(x). If e is negative, it defaults to
	//machine accuracy.
	Objective(double e, const Matrix& x0);
	//returns the current stored value of the objective funciton
	double getValue() const;
	//returns a copy of the current point of the objective, i.e. the argument
	//of the current value
	Matrix getArg() const;
	//returns the stored value of epsi. this is used for tolerance checking
	//in trust region and line search algorithms.
	double getEpsi() const;
	//Update the objective to a new point. Updates the value of the 
	//objective as well.
	void update(const Matrix& y);
	//An objective is defined as a scalar valued map. We use a matrix to
	//store the coordinates of the domain. The default objective is an
	//inner product.
	virtual double operator()(const Matrix& y);
	//We give a default numerical implementation of the gradient.
	//This just does an order 1 perturbation of every entry of x.
	//Any user of this library is free to overload the definition
	//and provide a different formula for gradient.
	//the result has the same shape as this->x.
	//If the gradient has already been computed for this x, then a flag
	//stating that the gradient is stored will be set to 1.
	//if that flag is 1, return the stored grad. Otherwise compute the
	//grad, store it, set flag to 1, and return.
	virtual Matrix& grad();
	//under construction:
		//We give a default numerical approximation of hessian
		//at a point x times a vector p.
		//This just calls the grad(x) function and computes an approximate
		//directional derivative via (grad(x + p*d) - grad(x))/d
		//virtual Matrix hessTimes(const Matrix& p);
	//This computes a full numerical approximation of the hessian
	//by finite differencing. A user of this library is free to overwrite
	//this definition and provide a different approximation or exact expression.
	//We don't store the hessian locally like the grad.
	virtual Matrix hess();
};

/*
 * Function Declarations:
 */
/*
 *	Does one iteration of the wolfe conditions line search algorithm to give
 *	a sufficient decrease update.
 * 
 */
Matrix wolfeLineSearch(	Objective& J,
							const Matrix& searchDirect,
							double c1,
							double c2,
							double alphaMax,
							int& stopSignal);

/*
 *	The zoom function for wolfe line search as defined by Nocedal and Wright
 *	in their text on Numerical Optimization
 */
double wolfeZoom(	Objective& J,
					const Matrix& searchDirect,
					double phiLo,
					double phiHi,
					double phiPrimeZero,
					double c1,
					double c2,
					double alphaLo,
					double alphaHi,
					int& stopSignal);

/*
 *	Does one iteration othe trust region method to give a sufficient decrease
 *	update. The quadratic subproblem is solved via the dogleg path.
 * 
 */
Matrix trustRegion(	Objective& J,
					const Matrix& gradAppx,
					const Matrix& hessAppx,
					double deltaHat,
					double delta_k,
					double eta);

/*
 *	Find the dogleg update for a trust region optimization framework.
 */
Matrix doglegStep(	const Matrix& gradAppx,
					const Matrix& hessAppx,
					double delta_k);

/*
 *	The BFGS quasi-newton optimization algorithm using wolfe conditions line
 *	search. Return value is number of iterations.
 *
 *	Stop Signal Manual:
 *	 Value  |         Meaning
 * ----------------------------------------------------------------------------
 *   -1     |  Hessian Not Invertible
 *    0     |  Successful Exit.
 *    1     |  Wolfe Line Search Failure, search direction gives no decrease.
 *    2     |  Wolfe Line Search Failure, search interval converges to 0 width.
 */
lapack_int bfgsLineSearch(	Objective& J,
							double c1,
							double c2,
							double alphaMax,
							int& stopSignal);

/*
 *	The limited memory version of  BFGS quasi-newton optimization algorithm
 *	using wolfe conditions line search.
 *	
 *	Stop Signal Manual:
 *	 Value  |         Meaning
 * ----------------------------------------------------------------------------
 *   -1     |  Hessian Not Invertible
 *    0     |  Successful Exit.
 *    1     |  Wolfe Line Search Failure, search direction gives no decrease.
 *    2     |  Wolfe Line Search Failure, search interval converges to 0 width.
 */
lapack_int bfgsLimitedMem(	Objective& J,
							double c1,
							double c2,
							double alphaMax,
							int& stopSignal,
							lapack_int histMaxSize = 0);

/*
 *	A quasi-newton method similar in flavor to BFGS which uses trust region
 *	to update iterations rather than wolfe conditions line search.
 */
Matrix sr1TrustRegion(	Objective& J,
						//std::valarray<Matrix>& hessArray,
						double delta_0,
						double eta,
						double r);

}
#endif
