#include "mtxLib.h"

typedef struct tUKFpar
{
    int stLen;//lenght of state vector
    double alfa;//Range:[10e-4 : 1].Smaller alfa leads to a tighter (closer) selection of sigma-points,
    double beta;//Contain information about the prior distribution (for Gaussian, beta = 2 is optimal).
    double k; //tertiary scaling parameter, usual value 0.
    double lambda;
    sMatrixType * pEta;

}tUKFpar;

typedef struct tUKFvar
{
    sMatrixType * pChi;//sigma points[stLen × (2*stLen+1)]
    sMatrixType * pPsi;//f(Chi[i])  
    sMatrixType * pPx;
    sMatrixType * pPy;
    double Ymean;

}tUKFvar;

typedef struct tUKF
{
    tUKFpar par; 
    tUKFvar var;
}tUKF;