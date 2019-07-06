#ifndef NUMPLES_ROOTFINDER_H
#define NUMPLES_ROOTFINDER_H

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
#include"numOpt.h"
#ifndef NULL
  #define NULL 0
#endif


namespace numples{
/*
 * Object and Struct Definitions:
 */
class Residual : public Objective{
private:
public:
  //A reference to the grad matrix stored in Objective class.
  Matrix& residMatrix;
	//Default constructor sets epsi to machine accuracy, val=0, x to zero.
	Residual() : Objective(), residMatrix(Objective::gradMatrix){
  };
	//A constructor which allows you to set epsi
	Residual(double e) : Objective(e), residMatrix(Objective::gradMatrix){
  };
	//A constructor which allows you to set initial guess yourself.
	//val is set to this->operator()(x).
	Residual(const Matrix& x0) : Objective(x0), 
                              residMatrix(Objective::gradMatrix){};
  //Need to define some constructors
  //returns the current stored value of the merit funciton
  double getValue() const{return Objective::getValue();}
  //returns a copy of the current point of the objective, i.e. the argument
  //of the current value
  Matrix getArg() const{return Objective::getArg();};
  //returns the stored value of epsi. this is used for tolerance checking
  //in trust region and line search algorithms.
  double getEpsi() const{return Objective::getEpsi();};
  //Update the objective to a new point. Updates the value of the 
  //objective as well.
  void update(const Matrix& y){Objective::update(y);};
  //This operator function is just a call to the merit() function. It is here
  //for use with trust region and line search algorithms.
  virtual double operator()(const Matrix& y);
  //A rootfinding algorithm tries to find points where r(x) = 0. r() is 
  //called the residual function. It is vector valued.
  virtual Matrix resid(const Matrix& y);
  //A merit function is a function we perform line search or truct region
  //algorithms on to find a small residual. default implementation is
  //half of residual norm squared.
  virtual double merit(const Matrix& y);
  //The grad function is different in that it returns value of residual
  //function, i.e. a vector/matrix from resid()
  virtual Matrix& grad();
  //hessian is not hessian, but instead the jacobian. it may not be 
  //symmetric. default implementation is finite-differences approximation.
  //virtual Matrix hess();
};
/*
 * Function Declarations:
 */

/*
 * An implementation of the Broyden Rootfinding algorithm which only requires
 * O(n*n) floating point operations per cycle.
 *
 *	Stop Signal Manual:
 *	 Value  |         Meaning
 * ----------------------------------------------------------------------------
 *   -1     |  Jacobian Not Invertible
 *    0     |  Successful Exit.
 *    1     |  Wolfe Line Search Failure, search direction gives no decrease.
 *    2     |  Wolfe Line Search Failure, search interval converges to 0 width.
 */
lapack_int broydenLineSearch(	Residual& J,
							                double c1,
							                double c2,
							                double alphaMax,
							                int& stopSignal);
}//end namespace numples
#endif
