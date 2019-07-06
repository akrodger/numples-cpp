/*
 * C Source file Template by Bram Rodgers.
 * Original Draft Dated: 25, Feb 2018
 */

/*
 * Macros and Includes go here: (Some common ones included)
 */
#include "rootFind.h"
#ifndef NULL
	#define NULL 0
#endif

namespace numples{
/*
 * Locally used helper functions:
 */

/*
 * Static Local Variables:
 */

/*
 * Function Implementations: eg:
 * Object Class::foo(Object in1, Object& in2);
 */

double Residual::operator()(const Matrix& y){
  return this->merit(y);
}

Matrix Residual::resid(const Matrix& y){
  return y;
}

double Residual::merit(const Matrix& y){
  Matrix r = this->resid(y);
  return 0.5*r.dot(r);
}

Matrix& Residual::grad(){
	if(this->gradAvailable != 0){
		return this->gradMatrix;
	}
  this->gradAvailable = 1;
  this->gradMatrix = this->resid(this->x);
  this->residMatrix = this->gradMatrix;
  return this->gradMatrix;
}

lapack_int broydenLineSearch(	Residual& J,
                              double c1,
                              double c2,
                              double alphaMax,
                              int& stopSignal){
  double iAngle = 0;
  lapack_int numIter= 0, probDim = 0;
  Matrix xNew;
  Matrix xOld = J.x;
  Matrix newResid;
  Matrix oldResid = J.grad();
  Matrix search;
  Matrix secant;
  Matrix residDif;
  Matrix zk = xOld;
  Matrix qk = xOld;
  Matrix appxJacInv;
  probDim = xOld.getNumComponents();//set j to dimension of x
  //Set initial Jacobian to be identity matrix.  
  appxJacInv = Matrix(probDim, probDim, 'I');
  while(1){
    search = appxJacInv*(oldResid*(-1.0));
    xNew = wolfeLineSearch( J,//Objective& J,
                            search,//const Matrix& searchDirect,
                            c1,//double c1,
                            c2,//double c2,
                            alphaMax,//double alphaMax,
                            stopSignal);//int& stopSignal);
    //Take new update and save to data structure
    J.update(xNew);
    newResid = J.grad();
    secant = xNew - xOld;
    residDif = newResid - oldResid;
    
    //check for halt condition
    if(newResid.norm2() < J.getEpsi() || (stopSignal != 0)){
      numIter++;
      break;
    }
    //Update the inverse jacobian with the scheme
    //H_new = H_old + u*v^T, with
    //update u,v vectors given by secant and residDif
    cblas_dgemv(CblasColMajor,
                CblasTrans,
                probDim,//OPENBLAS_CONST blasint m,
                probDim,//OPENBLAS_CONST blasint n,
                1.0,//OPENBLAS_CONST double alpha,
                &(appxJacInv.components[0]),//OPENBLAS_CONST double  *a,
                probDim,//OPENBLAS_CONST blasint lda,
                &(secant.components[0]),//OPENBLAS_CONST double  *x,
                1,//OPENBLAS_CONST blasint incx,
                0.0,//OPENBLAS_CONST double beta,
                &(zk.components[0]),//double  *y,
                1);//OPENBLAS_CONST blasint incy);
    iAngle = zk.dot(residDif);
    if(fabs(iAngle) < 1e-16){//check if inverse jacobian exists.
      stopSignal = -1;
      numIter++;
			break;
    }
    iAngle = 1.0 / iAngle;
    qk = appxJacInv*residDif;
    qk = secant - qk;
    cblas_dger(CblasColMajor,
                probDim,
                probDim,
                iAngle,
                &(zk.components[0]),
                1,
                &(qk.components[0]),
                1,
                &(appxJacInv.components[0]), probDim);
    //Done updating inverse jacobian.
    //update iterates and continue to next loop
    xOld = xNew;
    oldResid = newResid;
    numIter++;
  }
  return numIter;
}

};//end namespace numples
