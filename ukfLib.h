#include "mtxLib.h"

#define alphaIdx   (int)0
#define bethaIdx   (int)1
#define kappaIdx   (int)2
#define scalingLen  (int)3


typedef struct tUKFpar
{
    int stateLen;//lenght of state vector
    double alpha;//Range:[10e-4 : 1].Smaller alfa leads to a tighter (closer) selection of sigma-points,
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
    sMatrixType * pPk_sqrt; //1. pPcovXsqrt=sqrt(pPcovX)  &pPcovXsqrt = &pPcovX
    sMatrixType * pXk_p;
    sMatrixType * pSk_p;     // S(k-1)   Calculate the sigma-points
    sMatrixType * pSk_m;    // S(k|k-1) Propagate each sigma-point through prediction f(Chi)
    sMatrixType * pXk_m;    // X(k|k-1) Calculate mean of predicted state
    sMatrixType * pP_m;    // P(k|k-1) Calculate covariance of predicted state   
}tUKFpredict;

typedef struct tUKFobserv
{
    sMatrixType * pPsi_m; //1. Propagate each sigma-point through observation
    sMatrixType * pY_m;   //2. Calculate mean of predicted output
    sMatrixType * pPyy;   //3. Calculate covariance of predicted output
    sMatrixType * pPxy;   //4. Calculate cross-covariance of state and output
    
}tUKFobserv;

typedef struct tUKFupdate
{
    sMatrixType * pKgain;     //1. Calculate Kalman gain
    sMatrixType * pX;         //2. Update state estimate   Xk
    sMatrixType * pP;      //3. Update error covariance Pk


}tUKFupdate;

typedef struct tUKFvar
{
    tUKFpredict predict;
    tUKFobserv  observ;
    tUKFupdate  update;

    sMatrixType * pChi;//sigma points 
    sMatrixType * pPsi;//f(Chi[i])  

    sMatrixType * pPy;
    

}tUKFvar;

typedef struct tUKF
{
    tUKFpar par; 
    tUKFvar var;
}tUKF;