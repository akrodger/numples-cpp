/*
 * C Source file Template by Bram Rodgers.
 * Original Draft Dated: 25, Feb 2018
 */

/*
 * Macros and Includes go here: (Some common ones included)
 */
#include "numOpt.h"
#include<vector>
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

Objective::Objective(){
	//set epsi to machine accuracy
	this->epsi = 1.110223024625157e-16;
	//set x to zero
	this->x = Matrix(1,1,'0');
	//set val to zero
	this->val = 0;
	//set grad unavailable
	this->gradAvailable = 0;
	this->gradMatrix = Matrix(1,1,'0');
}

Objective::Objective(double e){
	//set epsi to machine accuracy
	this->epsi = e;
	//set x to zero
	this->x = Matrix(1,1,'0');
	//set val to zero
	this->val = 0;
	//set grad unavailable
	this->gradAvailable = 0;
	this->gradMatrix = Matrix(1,1,'0');
}

Objective::Objective(const Matrix& x0){
	//set epsi to machine accuracy
	this->epsi = 1.110223024625157e-16;
	//set x to input
	this->x = x0;
	//set val to objective of input
	this->val = this->operator()(x0);
	//set grad unavailable
	this->gradAvailable = 0;
	this->gradMatrix = Matrix(1,1,'0');
}

Objective::Objective(double e, const Matrix& x0){
	//set epsi to input
	this->epsi = e;
	//set x to input
	this->x = x0;
	//set val to input
	this->val = this->operator()(x0);
	//set grad unavailable
	this->gradAvailable = 0;
	this->gradMatrix = Matrix(1,1,'0');
}

double Objective::getValue() const{
	return this->val;
}

Matrix Objective::getArg() const{
	Matrix ecks = this->x;
	return ecks;
}

double Objective::getEpsi() const{
	return this->epsi;
}

void Objective::update(const Matrix& y){
	this->x = y;
	this->val = this->operator()(y);
	this->gradAvailable = 0;
}

double Objective::operator()(const Matrix& y){
	return y.dot(y);
}

Matrix Objective::grad(){
	lapack_int k = 0;
	lapack_int vecLen = this->x.getNumComponents();
	double perturb = 1.053671212772351e-08;
	//perturb is square root of epsi_mech, quasi optimal size
	double xHoldVal = 0;//holds unperturbed value of x in loop
	Matrix g;// = Matrix(this->x.getNumRows(), this->x.getNumCols());
	if(this->gradAvailable == 0){
		g = Matrix(this->x.getNumRows(), this->x.getNumCols());
		for(k = 0; k < vecLen; k++){
			xHoldVal = this->x(k);//hold onto value of x
			this->x(k) = this->x(k) + perturb;//perturb in direction k
			//get the value of the gradient in direction k
			g(k) = this->operator()(this->x);
			g(k) = (g(k) - this->getValue())/perturb;
			this->x(k) = xHoldVal;//return x to unperturbed value
		}
		this->gradAvailable = 1;
		this->gradMatrix = g;
	}else{
		return this->gradMatrix;
	}
	return g;
}

//Matrix Objective::hessTimes(const Matrix& p){
//	double perturb = 1.026484881901507e-04;//1.053671212772351e-08;
//	double yHoldVal = 0;
//	double yObjValue;
//	lapack_int k = 0;
//	lapack_int vecLen = this->x.getNumComponents();
//	Matrix y = this->x + (p*perturb);
//	Matrix Hp = Matrix(x.getNumRows(), x.getNumCols());
//	yObjValue = this->operator()(y);
//	for(k = 0; k < vecLen; k++){
//		yHoldVal = y(k);//hold onto value of y
//		y(k) = y(k) + perturb;//perturb in direction k
//		//get the value of the gradient in direction k
//		Hp(k) = (this->operator()(y) - yObjValue)/perturb;
//		y(k) = yHoldVal;//return x to unperturbed value
//	}
//	if(this->gradAvailable){
//		//finite diff to get directional hessian
//		Hp = (Hp - this->gradMatrix)/perturb;
//	}else{
//		Hp = (Hp - this->grad())/perturb;
//	}
//	return Hp;
//}

Matrix Objective::hess(){
	lapack_int j = 0, k = 0;
	lapack_int vecLen = this->x.getNumComponents();
	double perturb = 1.026484881901507e-04;//1.053671212772351e-08;
	//perturb is square root of epsi_mech, quasi optimal size
	double kHoldVal = 0;//holds unperturbed value of x in loop
	double jHoldVal = 0;//holds unperturbed value of x in loop
	Matrix H = Matrix(vecLen, vecLen);
	for(k = 0; k < vecLen; k++){
		kHoldVal = this->x(k);//hold onto value of x
		for(j = 0; j < k; j++){
			jHoldVal = this->x(j);//hold onto value of x
			this->x(k) = this->x(k) + perturb;//perturb in direction k
			this->x(j) = this->x(j) + perturb;//perturb in direction k
			//get the value of the gradient in direction k
			H(j,k) = this->operator()(this->x);
			this->x(k) = kHoldVal;//return x to unperturbed value
			H(j,k) = H(j,k) - this->operator()(this->x);
			this->x(j) = jHoldVal;
			this->x(k) = this->x(k) + perturb;
			H(j,k) = H(j,k) - this->operator()(this->x);
			this->x(k) = kHoldVal;//return x to unperturbed value
			H(j,k) = (H(j,k) + this->getValue())/(perturb*perturb);
			H(k,j) = H(j,k);
		}
		this->x(k) += perturb;
		H(k,k) = this->operator()(this->x);
		this->x(k) = kHoldVal - perturb;
		H(k,k) = H(k,k) +  this->operator()(this->x);
		H(k,k) = (H(k,k) - 2*this->getValue())/(perturb*perturb);
		this->x(k) = kHoldVal;
	}
	return H;
}

Matrix wolfeLineSearch(	Objective& J,
						const Matrix& searchDirect,
						double c1,
						double c2,
						double alphaMax,
						int& stopSignal){
	//Try full search length along searchdirect on first pass.
	double alpha_k = 1;//The full search length of searchDirect.
	double alpha_km1 = 0;
	double phiPrimeZero = searchDirect.dot(J.grad());
	double phiPrime_k = 0;
	double phi_k = 0;
	double phi_km1 = J.getValue();
	double alphaStar = 0;
	double perturb = 1.053671212772351e-08;
	int firstPass = -1;//first trial uses alpha_k = 1
	Matrix currPoint;
	while(1){
		//Do an update of objective along search direction
		currPoint = J.getArg() + (searchDirect*alpha_k);
		phi_k = J(currPoint);
		//If the sufficient decrease is satisfied or this iteration
		//gave an increase in the objective, then we know a strong
		//wolfe point exists. zoom and break.
		if((phi_k > (J.getValue() + c1*alpha_k*phiPrimeZero)) ||
			((phi_k >= phi_km1) && (firstPass == 0))){
			alphaStar = wolfeZoom(	J,
									searchDirect,
									phi_km1,
									phi_k,
									phiPrimeZero,
									c1,
									c2,
									alpha_km1,
									alpha_k,
									stopSignal);
			break;
		}
		//Do an update of the directional derivative of the objective
		phiPrime_k = J(J.getArg() + (searchDirect*(alpha_k + perturb)));
		phiPrime_k = (phiPrime_k - phi_k)/perturb;
		//if the curvature condition is satisfied, then then we can pick'
		//this iterate as our next optimizer point. so break.
		if(fabs(phiPrime_k) <= -c2*phiPrimeZero){
			alphaStar = alpha_k;
			break;
		}
		//If the derictional derivative points uphill, then we know that
		//there must have been a point were is was zero, by intermediate
		//value theorem. so zoom and then break.
		if(phiPrime_k >= 0){
			alphaStar = wolfeZoom(	J,
									searchDirect,
									phi_k,
									phi_km1,
									phiPrimeZero,
									c1,
									c2,
									alpha_k,
									alpha_km1,
									stopSignal);
			break;
		}
		//if alpha_k = 0.5*alphaMax didn't work,
		if(firstPass == 1){//then we are no longer in first trial case
			firstPass = 0;
		}
		//If we didn't get a full step on alpha_k = 1 trial,  
		if(firstPass == -1){//then we use a step between alpha and alphaMax
			alpha_k = alphaMax*0.5;
			alpha_km1 = 0;
			phi_k = 0;
			firstPass = 1;//now in original suggested first trial
			phi_km1 = J.getValue();
		}else{//we are in normal wolfe conditions deep iterations
			//We have yet to satisfy the wolfe conditions. cinch towards max
			//step size and iterate forward
			alpha_km1 = alpha_k;
			phi_km1 = phi_k;
			alpha_k = (alpha_k + alphaMax)*0.5;
		}
		//if step size is smaller than our differentiation perturbations,
		//then no optimizer update is reliable. halt.
		if(fabs(alpha_k - alpha_km1) < perturb){
			stopSignal = 1;
			break;
		}
		stopSignal = 0;//if nothing is wrong, then we can keep optimizing
	}
	currPoint = J.getArg() + (searchDirect*alphaStar);
	return currPoint;
}

double wolfeZoom(	Objective& J,
					const Matrix& searchDirect,
					double phiLo,
					double phiHi,
					double phiPrimeZero,
					double c1,
					double c2,
					double alphaLo,
					double alphaHi,
					int& stopSignal){
	double alpha_j = 0;
	double phi_j = 0;
	double phiPrime_j = 0;
	double perturb = 1.053671212772351e-08;
	double alphaStar = 0;
	Matrix currPoint;
	while(1){
		//Check if this iterval is big enough to be reliable for optimization
		if(fabs(alphaHi - alphaLo) > perturb){
			alpha_j = (alphaLo + alphaHi)*0.5;
		}else{//if it is not, then halt.
			stopSignal = 1;
			return (alphaLo + alphaHi)*0.5;
		}
		//Update objective along search direction.
		currPoint = J.getArg() + searchDirect*alpha_j;
		phi_j = J(currPoint);
		//By constriction, alphaLo satisfies sufficient decrease.
		//check to see if holds for phi_j or if phi_j is large enough
		//to ensure a directional minimizer between alphaLo and alpha_j
		if((phi_j > J.getValue() + c1*alpha_j*phiPrimeZero) ||
				(phi_j >= phiLo)){
			//if it does, then reel the high estimate backwards.
			alphaHi = alpha_j;
			phiHi = phi_j;
		}else{//If it does not, then modify points using curvature condition
			//first update the directional derivative
			phiPrime_j = J(J.getArg() + (searchDirect*(alpha_j + perturb)));
			phiPrime_j = (phiPrime_j - phi_j)/perturb;
			//check if the curvature condition holds. if it does, we're done
			if(fabs(phiPrime_j) <= -c2*phiPrimeZero){
				alphaStar = alpha_j;
				break;
			}
			//otherwise, double check that doesn't point uphill
			if(phiPrime_j*(alphaHi - alphaLo) >= 0){
				//if it does, reel alphaHi to the smaller value, alphaLo
				alphaHi = alphaLo;
				phiHi = phiLo;
			}
			//update the low value with our new alpha_j
			alphaLo = alpha_j;
			phiLo = phi_j;
		}
	}
	stopSignal = 0;
	return alphaStar;
}

lapack_int bfgsLineSearch(	Objective& J,
							double c1,
							double c2,
							double alphaMax,
							int& stopSignal){
	std::valarray<double> appxHessInv;
	double iAngle = 0;
	lapack_int i = 0, j = 0, numIter= 0, probDim = 0, setIndex = 0;
	Matrix xNew;
	Matrix xOld = J.getArg();
	Matrix newGrad;
	Matrix oldGrad = J.grad();
	Matrix search = Matrix(xOld.getNumRows(), xOld.getNumCols());
	Matrix secant;
	Matrix gradDif;
	probDim = xOld.getNumComponents();//set j to dimension of x
	//hessian symmetric, store only half of it.
	appxHessInv.resize(probDim*(probDim+1)/2);
	//Set initial hessian to be identity matrix.
	setIndex = 0;
	for(j = 0; j < probDim; j++){
		for(i = 0; i < j; i++){
			appxHessInv[setIndex] = 0;
			setIndex++;
		}
		appxHessInv[setIndex] = 1;
		setIndex++;
	}
	while(1){
		//search = oldGrad;//start computation of hessian inverse
		//Multiply the approximate inverse hessian into gradient and negate.
		cblas_dspmv(CblasColMajor, CblasUpper,
                 probDim, -1, &(appxHessInv[0]),
                 &(oldGrad.components[0]), 1,
                 0, &(search.components[0]), 1);

		xNew = wolfeLineSearch(	J,//Objective& J,
								search,//const Matrix& searchDirect,
								c1,//double c1,
								c2,//double c2,
								alphaMax,//double alphaMax,
								stopSignal);//int& stopSignal);
		//update the objective value
		J.update(xNew);
		newGrad = J.grad();
		//check for halt condition
		if(newGrad.norm2() < J.getEpsi() || (stopSignal != 0)){
			numIter++;
			break;
		}
		//get bfgs update parameters
		secant = xNew - xOld;//change in position this iterate
		gradDif = newGrad - oldGrad;//change in direction this iterate
		iAngle = 1/secant.dot(gradDif);//angle between changes
		//Use the xOld vector to temporarily hold curvature data
		cblas_dspmv(CblasColMajor, CblasUpper,
                 probDim, iAngle, &(appxHessInv[0]),
                 &(gradDif.components[0]), 1,
                 0, &(xOld.components[0]), 1);
		//xOld now stores hessInverse times gradDif
		//Now do a rank 2 update to the inverse hessian approx.
		cblas_dspr2(CblasColMajor,
					CblasUpper,
					probDim,
					-1,
					&(secant.components[0]),//OPENBLAS_CONST double *X, 
					1,//OPENBLAS_CONST blasint incX,
					&(xOld.components[0]),//OPENBLAS_CONST double *Y,
					1,//OPENBLAS_CONST blasint incY,
					&(appxHessInv[0]));
		//Now update the iAngle variable with an additional inner product
		iAngle = (gradDif.dot(xOld)*iAngle) + iAngle;
		cblas_dspr(CblasColMajor,//OPENBLAS_CONST enum CBLAS_ORDER order,
					CblasUpper,
					probDim,
					iAngle,//OPENBLAS_CONST double alpha,
					&(secant.components[0]),//OPENBLAS_CONST double *X,
					1,//OPENBLAS_CONST blasint incX,
					&(appxHessInv[0]));//double *Ap);
		//Done updating hessian.
		//update iterates and continue to next loop
		xOld = xNew;
		oldGrad = newGrad;
		numIter++;
	}//(newGrad.norm2() >= J.getEpsi() && (stopSignal == 0));
	return numIter;
}


lapack_int bfgsLimitedMem(	Objective& J,
							double c1,
							double c2,
							double alphaMax,
							int& stopSignal,
							lapack_int histMaxSize){
	std::vector<Matrix> secantHist;//history of secants (xNew - xOld)
	std::vector<Matrix> gradieHist;//history of gradient differences.
	std::vector<double> iAngleHist;//history of inverse search-to-grad angles
	std::vector<double> grProjHist;//history of secant projections onto search
	double lipschitz = 1;//An estimation of the lipshitz constant
	double betaProj = 1;//projection for second loop of hessian mult
	lapack_int histSize = 0, j = 0, numIter= 0;
	Matrix xNew;
	Matrix xOld = J.getArg();
	Matrix newGrad;
	Matrix oldGrad = J.grad();
	Matrix search;
	//If we are storing half of a matrix worth of history, then just
	//do the faster convering bfgs algorithm.
	if(histMaxSize >= xOld.getNumComponents()/4 || histMaxSize <= 0){
		return bfgsLineSearch(J, c1, c2, alphaMax, stopSignal);
	}
	while(1){
		search = oldGrad;//start computation of hessian inverse
		histSize = (lapack_int) secantHist.size();
		grProjHist.resize(histSize);
		//This loop computes the right half of the hessian inverse product
		for(j = histSize-1; j >= 0; j--){
			grProjHist[j] = (secantHist[j].dot(search))*iAngleHist[j];
			search = search - (gradieHist[j] * grProjHist[j]);
		}
		//this conditional does lipschitz constant multiplication.
		if(histSize > 0){
			j = histSize - 1;
			lipschitz = gradieHist[j].norm2();
			lipschitz *= lipschitz;
			lipschitz = secantHist[j].dot(gradieHist[j]);
			search = search*(lipschitz);
		}
		//this loop computes left half of inverse hessian multiplication
		for(j = 0; j < histSize; j++){
			betaProj = iAngleHist[j]*(gradieHist[j].dot(search));
			search = search + secantHist[j]*(grProjHist[j] - betaProj);
		}
		//point the search downhill
		search = search*-1;
		//wolfe conditions for next iterate
		xNew = wolfeLineSearch(	J,//Objective& J,
								search,//const Matrix& searchDirect,
								c1,//double c1,
								c2,//double c2,
								alphaMax,//double alphaMax,
								stopSignal);//int& stopSignal);
		//update the objective value
		J.update(xNew);
		newGrad = J.grad();
		//check for halt condition
		if(newGrad.norm2() < J.getEpsi() || (stopSignal != 0)){
			numIter++;
			break;
		}
		//load history of curvature data into storage vectors, pop old first
		if(histSize >= histMaxSize){
			secantHist.erase(secantHist.begin());
			gradieHist.erase(gradieHist.begin());
			iAngleHist.erase(iAngleHist.begin());
		}
		//save new curvature data
		secantHist.push_back(xNew - xOld);
		gradieHist.push_back(newGrad - oldGrad);
		j = ((lapack_int) secantHist.size()) - 1;
		iAngleHist.push_back(1 / (secantHist[j].dot(gradieHist[j])));
		//update iterates and continue to next loop
		xOld = xNew;
		oldGrad = newGrad;
		numIter++;
	}
	return numIter;
}

};
