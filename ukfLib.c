#include "ukfLib.h"
#include "stdio.h" //?ok or not ok
#include "math.h"


void ukf_init(tUKF * const pUkf,double scaling[scalingLen],int xLen,int yLen, tUkfMatrix * pUkfMatrix);


void ukf_step(tUKF * const pUkf);
void ukf_predict(tUKF * const pUkf);
void ukf_meas_update(tUKF * const pUkf);
void ukf_sigmapoint(tUKF * const pUkf);



void ukf_init(tUKF * const pUkf,double scaling[scalingLen],int xLen,int yLen, tUkfMatrix * pUkfMatrix)
{
    tUKFpar * pPar = (tUKFpar *)&pUkf->par;
    tUKFprev * pPrev = (tUKFprev *)&pUkf->prev;
    const int WmLen = pUkfMatrix->Wm_weight_vector->ncol;
    const int WcLen = pUkfMatrix->Wc_weight_vector->ncol;
    const int expWmLen = 2*xLen+1;
    
    mtx_cpy_f64(pUkf->prev.pP_p, pPar->pQ);
    //mtx_cpy_f64(pUkf->prev.pP_p, pPar->pR);

    pPar->pWm = pUkfMatrix->Wm_weight_vector;
    pPar->alpha = scaling[alphaIdx];
    pPar->betha = scaling[bethaIdx];
    pPar->kappa = scaling[kappaIdx];
    pPar->xLen = xLen;
    pPar->yLen = yLen;//check length of output matrix to compare
    pPar->sLen = 2*xLen+1;
    

    //calculate UKF lambda parameter
    pPar->lambda = pPar->alpha * pPar->alpha;
    pPar->lambda *= (double)(xLen + pPar->kappa);
    pPar->lambda -= (double)xLen;
    
    if(WmLen == expWmLen && WcLen == WmLen)
    {
        int col;
        const double Wm0 = pPar->lambda/(pPar->xLen + pPar->lambda);

        pPar->pWm->val[0] = Wm0;
        pPar->pWc->val[0] = Wm0 + (1 - pPar->alpha*pPar->alpha + pPar->betha) ;

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
    pPrev->pX_p = pUkfMatrix->X_sigma_points; //share same memory with pX_m
    pPrev->pu_p = pUkfMatrix->u_system_input;
    pPrev->px_p = pUkfMatrix->x_system_states;

    pUkf->predict.pP_m = pUkfMatrix->P_error_covariance;
    pUkf->predict.pX_m = pUkfMatrix->X_sigma_points;
    pUkf->predict.px_m = pUkfMatrix->x_system_states;
    pUkf->predict.pY_m = pUkfMatrix->Y_sigma_points;
    pUkf->predict.py_m = pUkfMatrix->y_predicted_mean;

    pUkf->update.pIxx = pUkfMatrix->I_identity_matrix;
    pUkf->update.pK = pUkfMatrix->K_kalman_gain;
    pUkf->update.pP = pUkfMatrix->P_error_covariance;
    pUkf->update.pPxy = pUkfMatrix->Pxy_cross_covariance;
    pUkf->update.pPyy = pUkfMatrix->Pyy_out_covariance;
    pUkf->update.px = pUkfMatrix->x_system_states; //&px = &px_m = &px_p
}

void ukf_step(tUKF * const pUkf)
{
    ukf_sigmapoint(pUkf);
    ukf_predict(pUkf);
    ukf_meas_update(pUkf);  
}

//Step 1: Generate the Sigma-Points
void ukf_sigmapoint(tUKF * const pUkf)
{
    double * pP_p = (double *)&pUkf->prev.pP_p->val;
    double * pX_p = (double *)&pUkf->prev.pX_p->val;
    double * px_p = (double *)&pUkf->prev.px_p->val;
    double scalarMultiplier;
    const double lambda = pUkf->par.lambda;
    const int sLen = pUkf->par.sLen;
    const int xLen = pUkf->par.xLen;
    int xIdx;
    int sigmaIdx=0;
    
    //1. Calculate error covariance matrix square root: P_p = sqrt(P_p)
    (void)mtx_chol_f64(pUkf->prev.pP_p);
    
    //2. Calculate the sigma-points
    for(xIdx=0;xIdx<xLen;xIdx++)
    {
        //first column of sigma point matrix is equal of previous state value
        pX_p[sLen*xIdx+sigmaIdx] = px_p[xIdx];
    }
    
    scalarMultiplier = sqrt(xLen + lambda);
    
    (void)mtx_mul_scalar_f64(pUkf->prev.pP_p,scalarMultiplier);
    
    for(sigmaIdx=1;sigmaIdx < sLen;sigmaIdx++)
    {
        for(xIdx=0;xIdx<xLen;xIdx++)
        {
            if(sigmaIdx <= xLen)
            {
                pX_p[sLen*xIdx+sigmaIdx] = px_p[xIdx] + pP_p[xLen*xIdx + (sigmaIdx-1)];
            }
            else
            {
                pX_p[sLen*xIdx+sigmaIdx] = px_p[xIdx] - pP_p[xLen*xIdx + (sigmaIdx-5)];
            }            
        }
    }    
}

//Step 2: Prediction Transformation
void ukf_predict(tUKF * const pUkf)
{
    tPredictFcn * pFcnPredict = (tPredictFcn *)&pUkf->predict.pFcnPredict;
    tObservFcn * pFcnObserve = (tObservFcn *)&pUkf->predict.pFcnObserv;
    tUKFpar * pPar = (tUKFpar *)&pUkf->par;
    double * pWm = (double *)&pPar->pWm->val;
    double * pWc = (double *)&pPar->pWc->val;
    double * pX_m = (double *)&pUkf->predict.pX_m->val;
    double * pY_m = (double *)&pUkf->predict.pY_m->val;
    double * pP_m = (double *)&pUkf->predict.pP_m->val;
    double * pPyy = (double *)&pUkf->update.pPyy;
    double * pPxy = (double *)&pUkf->update.pPxy;
    double * px_m = (double *)&pUkf->predict.px_m;
    double * py_m = (double *)&pUkf->predict.py_m;
    const int sigmaLen = pPar->sLen;
    const int xLen = pPar->xLen;
    const int yLen= pPar->yLen;
    int sigmaIdx,xIdx,xTrIdx,yIdx,yTrIdx;
    

    for(xIdx=0;xIdx<xLen;xIdx++)
    {
        px_m[xIdx] = 0;
               
        for(sigmaIdx=0;sigmaIdx<sigmaLen;sigmaIdx++)
        {                     
            if(pFcnPredict[xIdx] != NULL)
            {
                //Propagate each sigma-point through prediction
                //X_m[L][2L+1] = f(X_p, u_p) 
                pFcnPredict[xIdx](pUkf->prev.pu_p, pUkf->prev.pX_p, pUkf->predict.pX_m,sigmaIdx);
            }

            //Calculate mean of predicted state 
            //x_m[L][1] = sum(Wm(i)*X_m(i))
            px_m[xIdx] += pWm[sigmaIdx] * pX_m[sigmaLen*sigmaIdx + xIdx];
        }
    }
    
    //loop for each sigma point column vector  
    for(sigmaIdx=0;sigmaIdx<sigmaLen;sigmaIdx++)
    {       
        //3.Calculate covariance of predicted state P(k|k-1)
        //loop row of COV[:][L]
        for(xIdx=0;xIdx<xLen;xIdx++)
        {                      
            //ToDo: P(k|k-1) = Q(k-1) 
            for(xTrIdx=0;xTrIdx<xLen;xTrIdx++)
            {
                //loop col of COV[L][:] 
                
                double term1 = (pX_m[sigmaLen*xIdx + sigmaIdx] - px_m[xIdx]);
                double term2 = (pX_m[sigmaLen*xTrIdx + sigmaIdx] - px_m[xIdx]);
                
                //Perform multiplication with accumulation for each covariance matrix index
                //P(k|k-1)[xIdx][xTrIdx]= Wc(sigmaIdx)*(X_m-x_mean)*(X_m-x_mean)'
                //result is 2L+1 COV[L][L] matrix which are weighted in one common
                pP_m[xLen*xIdx + xTrIdx] += pWc[sigmaIdx]*term1*term2;
            }           
        } 
        
        //Calculate outputs for each sigma point
        py_m[yIdx] = 0;
        for(yIdx=0;yIdx<yLen;yIdx++)
        {
            if(pFcnObserve[yIdx] != NULL)
            {
                //Propagate each sigma-point through observation Y_m[yL][2L+1] = h(X_m, u)             
                pFcnObserve[yIdx](pUkf->input.pu, pUkf->predict.pX_m,pUkf->predict.pY_m,sigmaIdx);
            }
            else
            {
                //assign 0 if observation function is not specified
                pY_m[sigmaLen*sigmaIdx + yIdx]=0;
            }
            //Calculate mean of predicted output 
            //y_m[L] = sum(Wm(i)*Y_m(i))
            py_m[yIdx] += pWm[sigmaIdx] * pY_m[sigmaLen*sigmaIdx + yIdx];
        }
    }

    for(sigmaIdx=0;sigmaIdx<sigmaLen;sigmaIdx++)
    {       
        //Calculate covariance of predicted output
        //loop row of Pyy[:][yL]
        for(yIdx=0;yIdx<yLen;yIdx++)
        {                      
            //ToDo: Pyy(k|k-1) = R(k) 
            for(yTrIdx=0;yTrIdx<yLen;yTrIdx++)
            {
                //loop col of COV[L][:] 
                
                double term1 = (pY_m[sigmaLen*yIdx + sigmaIdx] - py_m[yIdx]);
                double term2 = (pY_m[sigmaLen*yTrIdx + sigmaIdx] - py_m[yIdx]);
                
                //Perform multiplication with accumulation for each covariance matrix index
                //Pyy(k)[yIdx][yTrIdx]= Wc(sigmaIdx)*(Y_m-y_m)*(Y_m-y_m)'
                //result is (2L+1)x(Pyy[yL][yL]) matrix which are weighted in one common
                pPyy[yLen*yIdx + yTrIdx] += pWc[sigmaIdx]*term1*term2;            
            }           
        }
        
        //Calculate cross-covariance of state and output
        for(xIdx=0;xIdx<xLen;xIdx++)
        {
            for(yTrIdx=0;yTrIdx<yLen;yTrIdx++)
            {
                double term1 = (pX_m[sigmaLen*xIdx + sigmaIdx] - px_m[xIdx]);
                double term2 = (pY_m[sigmaLen*yTrIdx + sigmaIdx] - py_m[yIdx]);

                
                pPxy[yLen*xIdx + yTrIdx] += pWc[sigmaIdx]*term1*term2;
            }
        }
    }
}

void ukf_meas_update(tUKF * const pUkf)
{
    tUKFupdate * pUpdate = (tUKFupdate*)&pUkf->update;


    //3.Calculate Kalman gain: 
    (void)mtx_identity_f64(pUpdate->pIxx);
    
    //inv(Pyy)
    (void)mtx_inv_f64(pUpdate->pPyy,pUpdate->pIxx);

    //K = Pxy * inv(Pyy)
    (void)mtx_mul_f64(pUpdate->pPxy,pUpdate->pIxx,pUpdate->pK);
    
    //K'
    (void)mtx_transp_dest_f64(pUpdate->pK, pUpdate->pKt);

    //4.Update state estimate

    // y = y - y_m
    (void)mtx_subtract_f64(pUkf->input.py,pUkf->predict.py_m);

    // K*(y - y_m) temporal stored in py
    (void)mtx_mul_f64(pUpdate->pK, pUkf->input.py, pUkf->input.py);

    // x_m + K*(y - y_m)
    (void)mtx_add_f64(pUkf->predict.px_m, pUkf->input.py);


    //5.Update error covariance
    //use Pxy for temporal resul trom multiplication
    //Pxy = K*Pyy
    (void)mtx_mul_f64(pUpdate->pK,pUpdate->pPyy,pUpdate->pPxy);

    //Ixx = K*Pyy*K'
    (void)mtx_mul_f64(pUpdate->pPxy, pUpdate->pKt, pUpdate->pIxx);

    (void)mtx_subtract_f64(pUkf->predict.pX_m,pUpdate->pIxx);

}




