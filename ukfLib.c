#include "ukfLib.h"
#include "stdio.h" //?ok or not ok


void ukf_init(tUKF * const pUkf,double scaling[scalingLen],int xLen, tUkfMatrix * pUkfMatrix);


void ukf_step(tUKF * const pUkf);
void ukf_predict(tUKF * const pUkf);
void ukf_observe(tUKF * const pUkf);
void ukf_meas_update(tUKF * const pUkf);
void ukf_sigmapoint(tUKF * const pUkf);



void ukf_init(tUKF * const pUkf,double scaling[scalingLen],int xLen, tUkfMatrix * pUkfMatrix)
{
    tUKFpar * pPar = (tUKFpar *)&pUkf->par;
    tUKFprev * pPrev = (tUKFprev *)&pUkf->prev;
    const int WmLen = pUkfMatrix->Wm_weight_vector->ncol;
    const int WcLen = pUkfMatrix->Wc_weight_vector->ncol;
    const int expWmLen = 2*xLen+1;
    

    pPar->pWm = pUkfMatrix->Wm_weight_vector;
    pPar->alpha = scaling[alphaIdx];
    pPar->betha = scaling[bethaIdx];
    pPar->kappa = scaling[kappaIdx];
    pPar->xLen = xLen;
    

    //calculate UKF lambda parameter
    pPar->lambda = pPar->alpha * pPar->alpha;
    pPar->lambda *= (double)(xLen + pPar->kappa);
    pPar->lambda -= (double)xLen;
    
    if(WmLen == expWmLen && WcLen == WmLen)
    {
        int col = 0;

        pPar->pWm->val[col] = 0;
        pPar->pWc->val[col] = 0;

        for(col=1;col<WmLen;col++)
        {
            pPar->pWm->val[col] = 1 / (2*(pPar->xLen + pPar->lambda));
            pPar->pWc->val[col] = pPar->pWm->val[col];
        }
    }
    else
    {
        //UKF init fail 
    }

    pUkf->input.pu = pUkfMatrix->u_system_input;

    pPrev->pP_p = pUkfMatrix->P_error_covariance;
    pPrev->pX_p = pUkfMatrix->X_sigma_point_gen; //share same memory with pX_m
    pPrev->pu_p = pUkfMatrix->u_system_input;
    pPrev->px_p = pUkfMatrix->x_system_states;

    pUkf->predict.pP_m = pUkfMatrix->P_error_covariance;
    pUkf->predict.pX_m = pUkfMatrix->X_sigma_point_gen;
    pUkf->predict.px_m = pUkfMatrix->x_system_states;
    pUkf->predict.pY_m = pUkfMatrix->Y_sigma_point_propagated;
    pUkf->predict.py_m = pUkfMatrix->y_predicted_mean;

    pUkf->update.pK = pUkfMatrix->K_kalman_gain;
    pUkf->update.pP = pUkfMatrix->P_error_covariance;
    pUkf->update.pPxy = pUkfMatrix->Pxy_cross_covariance;
    pUkf->update.pPyy = pUkfMatrix->Pyy_predicted_covariance;
    pUkf->update.px = pUkfMatrix->x_system_states; //&px = &px_m = &px_p
}

void ukf_step(tUKF * const pUkf)
{
    ukf_sigmapoint(pUkf);
    ukf_predict(pUkf);
    ukf_observe(pUkf);
    ukf_meas_update(pUkf);  
}

//Step 1: Generate the Sigma-Points
void ukf_sigmapoint(tUKF * const pUkf)
{
    //1. Calculate error covariance matrix square root

    //2. Calculate the sigma-points
    
}

//Step 2: Prediction Transformation
void ukf_predict(tUKF * const pUkf)
{
    tPredictFcn * pFcnPredict = (tPredictFcn *)&pUkf->predict.pFcnPredict;
    tUKFpar * pPar = (tUKFpar *)&pUkf->par;
    double * pWm = (double *)&pPar->pWm->val;
    double * pWc = (double *)&pPar->pWc->val;
    double * pX_m = (double *)&pUkf->predict.pX_m->val;
    double * px_m = (double *)&pUkf->predict.px_m;
    const int sigmaLen = 2*pPar->xLen+1;
    const int stateLen = pPar->xLen;
    int sigmaIdx,stateIdx;
    

    for(stateIdx=0;stateIdx<stateLen;stateIdx++)
    {
        px_m[stateIdx] = 0;
        
        for(sigmaIdx=0;sigmaIdx<sigmaLen;sigmaIdx++)
        {                     
            if(pFcnPredict[stateIdx] != NULL)
            {
                //1. Propagate each sigma-point through prediction X_m[L][2L+1] = f(X_p, u_p) 
                pFcnPredict[stateIdx](pUkf->prev.pu_p, pUkf->prev.pX_p, pUkf->predict.pX_m,sigmaIdx);
            }
            //2.Calculate mean of predicted state x_m[L][1] = sum(Wm(i)*X_m(i))
            px_m[stateIdx] += pWm[sigmaIdx] * pX_m[sigmaLen*sigmaIdx + stateIdx];

        }
    }

    

    //3.Calculate covariance of predicted state P(k|k-1)
    //consider to upgarde lib

    for(sigmaIdx=0;sigmaIdx<sigmaLen;sigmaIdx++)
    {
        double * pP_m = (double *)&pUkf->predict.pP_m->val;
        
        for(stateIdx=0;stateIdx<stateLen;stateIdx++)
        {
            int stateTrIdx;

            for(stateTrIdx=0;stateTrIdx<stateLen;stateTrIdx++)
            {
                double term1 = (pX_m[sigmaLen*stateIdx + sigmaIdx] - px_m[stateIdx]);
                double term2 = (pX_m[sigmaLen*stateTrIdx + sigmaIdx] - px_m[stateIdx]);

                pP_m[stateLen*stateIdx + stateTrIdx] = pWm[sigmaIdx]*term1*term2;
            }           
        }
        //pP_m[L][L]
        //pP_m[sigmaLen*sigmaIdx + stateIdx] += pWm[sigmaIdx];

        
    }
    //4.Propagate each sigma-point through observation

    //5.Calculate mean of predicted output
}

void ukf_observe(tUKF * const pUkf)
{
    tUKFpar * pPar = (tUKFpar *)&pUkf->par;
    tUKFpredict * pPredict = (tUKFpredict *)&pUkf->predict;
    const int sigmaLen = 2*pPar->xLen+1;
    int sigmaIdx;
    
    
    for(sigmaIdx=0;sigmaIdx<sigmaLen;sigmaIdx++)
    {
        //1.Propagate each sigma-point through observation
        pPredict->pFcnObserv[sigmaIdx](pUkf->input.pu, pUkf->predict.pX_m, pUkf->predict.pY_m,sigmaIdx);

        //2.Calculate mean of predicted output
    }



}

void ukf_meas_update(tUKF * const pUkf)
{
    tUKFupdate * pUpdate = (tUKFupdate*)&pUkf->update;

    //1.Calculate covariance of predicted output

    //2.Calculate cross-covariance of state and output

    //3.Calculate Kalman gain: Pxy * inv(Pyy)

    (void)mtx_inv_f64(pUpdate->pPyy,pUpdate->pI);
    (void)mtx_mul_f64(pUpdate->pPxy,pUpdate->pI,pUpdate->pK);
    
    //4.Update state estimate

    //5.Update error covariance

}




