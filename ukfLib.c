#include "ukfLib.h"


void ukf_init(tUKF * const pUkf,double alfa,double beta,double k,int stLen, sMatrixType * pCovX, sMatrixType * pCovY, sMatrixType * pSigma, sMatrixType * pEta);


void ukf_init(tUKF * const pUkf,double alfa,double beta,double k,int stLen, sMatrixType * pCovX, sMatrixType * pCovY, sMatrixType * pSigma, sMatrixType * pEta)
{
    tUKFpar * ukfPar = (tUKFpar *)&pUkf->par;
    tUKFvar * ukfVar = (tUKFvar *)&pUkf->var;
    const int nRows = pEta->nrow;
    int row;

    ukfPar->pEta = pEta;
    ukfPar->alfa = alfa;
    ukfPar->beta = beta;
    ukfPar->stLen = stLen;
    ukfPar->k = k;

    ukfPar->lambda = alfa * alfa;
    ukfPar->lambda *= (double)(stLen + k);
    ukfPar->lambda -= (double)stLen;

    for(row=0;row<nRows/* 2*stLen */;row++)
    {
        ukfPar->pEta->val[row] = 1 / (2*(stLen + ukfPar->lambda));
    }


    ukfVar->pPx = pCovX;
    ukfVar->pPy = pCovY;
    ukfVar->pChi = pSigma;

}


