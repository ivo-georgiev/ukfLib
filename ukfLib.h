#include "mtxLib.h"

#define alphaIdx   (int)0
#define bethaIdx   (int)1
#define kappaIdx   (int)2
#define scalingLen  (int)3


typedef struct tUKFpar
{
    int stateLen;//lenght of state vector
    double alpha;//Range:[10e-4 : 1].Smaller alpha leads to a tighter (closer) selection of sigma-points,
    double betha;//Contain information about the prior distribution (for Gaussian, beta = 2 is optimal).
    double kappa; //tertiary scaling parameter, usual value 0.
    double lambda;
    sMatrixType * pWm;
    sMatrixType * pWc;
    sMatrixType * pQ;
    sMatrixType * pR;

}tUKFpar;

typedef struct tUKFpredict //p(previous)==k-1, m(minus)=(k|k-1)
{
    sMatrixType * pPsqrt; //1. pPcovXsqrt=sqrt(pPcovX)  &pPcovXsqrt = &pPcovX
    sMatrixType * px_p;    // x(k-1)
    sMatrixType * pX_p;    // X(k-1)   Calculate the sigma-points
    sMatrixType * pX_m;    // X(k|k-1) Propagate each sigma-point through prediction f(Chi)
    sMatrixType * px_m;    // x(k|k-1) Calculate mean of predicted state
    sMatrixType * pP_m;    // P(k|k-1) Calculate covariance of predicted state   
}tUKFpredict;

typedef struct tUKFobserv
{
    sMatrixType * pY_m;   //1. Y(k|k-1) Propagate each sigma-point through observation
    sMatrixType * py_m;   //2. y(k|k-1) Calculate mean of predicted output
    sMatrixType * pPyy;   //3. Calculate covariance of predicted output
    sMatrixType * pPxy;   //4. Calculate cross-covariance of state and output
    
}tUKFobserv;

typedef struct tUKFupdate
{
    sMatrixType * pK;     //1. K(k) Calculate Kalman gain
    sMatrixType * px;     //2. x(k) Update state estimate   
    sMatrixType * pP;     //3. P(k) Update error covariance 


}tUKFupdate;

typedef struct tUKFvar
{
    tUKFpredict predict;
    tUKFobserv  observ;
    tUKFupdate  update;
}tUKFvar;

typedef struct tUKF
{
    tUKFpar par; 
    tUKFvar var;
}tUKF;