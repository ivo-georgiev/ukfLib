#include "ukfLib.h"


void ukf_init(tUKF * const pUkf,double scaling[scalingLen],int stateLen, sMatrixType * pCovX, sMatrixType * pCovY, sMatrixType * pSigma, sMatrixType * pWm, sMatrixType * pWc);


void ukf_init(tUKF * const pUkf,double scaling[scalingLen],int stateLen, sMatrixType * pCovX, sMatrixType * pCovY, sMatrixType * pSigma, sMatrixType * pWm, sMatrixType * pWc)
{
    tUKFpar * ukfPar = (tUKFpar *)&pUkf->par;
    tUKFvar * ukfVar = (tUKFvar *)&pUkf->var;
    const int WmLen = pWm->ncol;
    const int WcLen = pWc->ncol;
    const int expWmLen = 2*stateLen+1;
    

    ukfPar->pWm = pWm;
    ukfPar->alpha = scaling[alphaIdx];
    ukfPar->betha = scaling[bethaIdx];
    ukfPar->kappa = scaling[kappaIdx];
    ukfPar->stateLen = stateLen;
    

    //calculate UKF lambda parameter
    ukfPar->lambda = ukfPar->alpha * ukfPar->alpha;
    ukfPar->lambda *= (double)(stateLen + ukfPar->kappa);
    ukfPar->lambda -= (double)stateLen;
    
    if(WmLen == expWmLen && WcLen == WmLen)
    {
        int col = 0;

        ukfPar->pWm->val[col] = 0;
        ukfPar->pWc->val[col] = 0;

        for(col=1;col<WmLen;col++)
        {
            ukfPar->pWm->val[col] = 1 / (2*(ukfPar->stateLen + ukfPar->lambda));
            ukfPar->pWc->val[col] = ukfPar->pWm->val[col];
        }
    }
    else
    {
        //UKF init fail 
    }



//    ukfVar->predict.pPk= pCovX;
//    ukfVar->pPy = pCovY;
//    ukfVar->pChi = pSigma;

}


