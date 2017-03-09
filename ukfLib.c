#include "ukfLib.h"
#include "stdio.h" //?ok or not ok
#include "math.h"
#include "memory.h"


void ukf_init(tUKF * const pUkf,double scaling[scalingLen],int xLen,int yLen, tUkfMatrix * pUkfMatrix);


void ukf_step(tUKF * const pUkf);
void ukf_predict(tUKF * const pUkf);
void ukf_meas_update(tUKF * const pUkf);
void ukf_sigmapoint(tUKF * const pUkf);



void ukf_init(tUKF * const pUkf,double scaling[scalingLen],int xLen,int yLen, tUkfMatrix * pUkfMatrix)
{
    tUKFpar * pPar = (tUKFpar *)&pUkf->par;
    tUKFprev * pPrev = (tUKFprev *)&pUkf->prev;
    const int WmLen = pUkfMatrix->Wm_weight_vector.ncol;
    const int WcLen = pUkfMatrix->Wc_weight_vector.ncol;
    const int expWmLen = 2*xLen+1;
    
    pPar->Pxx0 = pUkfMatrix->Pxx0_init_error_covariance;
    pPar->Qxx = pUkfMatrix->Qxx_process_noise_cov;
    pPar->Wm =  pUkfMatrix->Wm_weight_vector;
    pPar->Wc =  pUkfMatrix->Wc_weight_vector;
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

        pPar->Wm.val[0] = Wm0;
        pPar->Wc.val[0] = Wm0 + (1 - pPar->alpha*pPar->alpha + pPar->betha) ;

        for(col=1;col<WmLen;col++)
        {
            pPar->Wm.val[col] = 1 / (2*(pPar->xLen + pPar->lambda));
            pPar->Wc.val[col] = pPar->Wm.val[col];
        }
    }
    else
    {
        //UKF init fail 
    }

    pUkf->input.u = pUkfMatrix->u_system_input;

    pPrev->Pxx_p = pUkfMatrix->Pxx_error_covariance;
    mtx_cpy_f64(&pUkf->prev.Pxx_p, &pPar->Pxx0);
    //mtx_cpy_f64(pUkf->prev.pP_p, pPar->pR);

    pPrev->X_p = pUkfMatrix->X_sigma_points; //share same memory with X_m
    pPrev->u_p = pUkfMatrix->u_system_input;
    pPrev->x_p = pUkfMatrix->x_system_states;

    pUkf->predict.P_m = pUkfMatrix->Pxx_error_covariance;
    pUkf->predict.X_m = pUkfMatrix->X_sigma_points;
    pUkf->predict.x_m = pUkfMatrix->x_system_states;
    pUkf->predict.Y_m = pUkfMatrix->Y_sigma_points;
    pUkf->predict.y_m = pUkfMatrix->y_predicted_mean;
    pUkf->predict.pFcnPredict = pUkfMatrix->fcnPredict;
    pUkf->predict.pFcnObserv = pUkfMatrix->fcnObserve;

    pUkf->update.Ixx = pUkfMatrix->I_identity_matrix;
    pUkf->update.K = pUkfMatrix->K_kalman_gain;
    pUkf->update.Pxx = pUkfMatrix->Pxx_error_covariance;
    pUkf->update.Pxy = pUkfMatrix->Pxy_cross_covariance;
    pUkf->update.Pyy = pUkfMatrix->Pyy_out_covariance;
    pUkf->update.x = pUkfMatrix->x_system_states; //&px = &px_m = &px_p

   
}

void ukf_step(tUKF * const pUkf)
{
    ukf_predict(pUkf);
    ukf_meas_update(pUkf);  
}

//Step 1: Generate the Sigma-Points
void ukf_sigmapoint(tUKF * const pUkf)
{
    double * Pxx_p = pUkf->prev.Pxx_p.val;
    double * X_p = pUkf->prev.X_p.val;
    double * x_p = pUkf->prev.x_p.val;
    double scalarMultiplier;
    const double lambda = pUkf->par.lambda;
    const int sLen = pUkf->par.sLen;
    const int xLen = pUkf->par.xLen;
    int xIdx;
    int sigmaIdx=0;
    mtxResultInfo mtxResult;
    
    //1. Calculate error covariance matrix square root: P_p = sqrt(P_p)
    mtxResult = mtx_chol_f64(&pUkf->prev.Pxx_p);

    if(MTX_OPERATION_OK == mtxResult)
    {       
        //2. Calculate the sigma-points
        for(xIdx=0;xIdx<xLen;xIdx++)
        {
            //first column of sigma point matrix is equal of previous state value
            X_p[sLen*xIdx+sigmaIdx] = x_p[xIdx];
        }
        
        scalarMultiplier = sqrt(xLen + lambda);
        
        (void)mtx_mul_scalar_f64(&pUkf->prev.Pxx_p,scalarMultiplier);
        
        for(sigmaIdx=1;sigmaIdx < sLen;sigmaIdx++)
        {
            for(xIdx=0;xIdx<xLen;xIdx++)
            {
                if(sigmaIdx <= xLen)
                {
                    X_p[sLen*xIdx+sigmaIdx] = x_p[xIdx] + Pxx_p[xLen*xIdx + (sigmaIdx-1)];//bug bug
                }
                else
                {
                    X_p[sLen*xIdx+sigmaIdx] = x_p[xIdx] - Pxx_p[xLen*xIdx + (sigmaIdx-5)];
                }            
            }
        } 
    }
    else
    {

    }
}

//Step 2: Prediction Transformation
void ukf_predict(tUKF * const pUkf)
{
    tUKFpar * pPar = (tUKFpar *)&pUkf->par;
    double * Wm = pPar->Wm.val;
    double * Wc = pPar->Wc.val;
    double * X_m = pUkf->predict.X_m.val;
    double * pY_m = pUkf->predict.Y_m.val;
    double * pP_m = pUkf->predict.P_m.val;
    double * pPyy = (double *)&pUkf->update.Pyy;
    double * Pxy = (double *)&pUkf->update.Pxy;
    double * px_m = pUkf->predict.x_m.val;
    double * py_m = pUkf->predict.y_m.val;
    const int sigmaLen = pPar->sLen;
    const int xLen = pPar->xLen;
    const int yLen= pPar->yLen;
    int sigmaIdx,xIdx,xTrIdx,yIdx,yTrIdx;
    
    ukf_sigmapoint(pUkf);

    for(xIdx=0;xIdx<xLen;xIdx++)
    {
        px_m[xIdx] = 0;
               
        for(sigmaIdx=0;sigmaIdx<sigmaLen;sigmaIdx++)
        {                     
            if(pUkf->predict.pFcnPredict[xIdx] != NULL)
            {
                //Propagate each sigma-point through prediction
                //X_m[L][2L+1] = f(X_p, u_p) 
                pUkf->predict.pFcnPredict[xIdx](&pUkf->prev.u_p, &pUkf->prev.X_p, &pUkf->predict.X_m,sigmaIdx);
            }

            //Calculate mean of predicted state 
            //x_m[L][1] = sum(Wm(i)*X_m(i))
            px_m[xIdx] += Wm[sigmaIdx] * X_m[sigmaLen*sigmaIdx + xIdx];
        }
    }

    //py_m[yIdx] = 0;
    memset(py_m,0,sizeof(double)*sigmaLen*yLen);
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
                
                double term1 = (X_m[sigmaLen*xIdx + sigmaIdx] - px_m[xIdx]);
                double term2 = (X_m[sigmaLen*xTrIdx + sigmaIdx] - px_m[xIdx]);
                
                //Perform multiplication with accumulation for each covariance matrix index
                //P(k|k-1)[xIdx][xTrIdx]= Wc(sigmaIdx)*(X_m-x_mean)*(X_m-x_mean)'
                //result is 2L+1 COV[L][L] matrix which are weighted in one common
                pP_m[xLen*xIdx + xTrIdx] += Wc[sigmaIdx]*term1*term2;
            }           
        } 
        
        //Calculate outputs for each sigma point
        
        for(yIdx=0;yIdx<yLen;yIdx++)
        {
            if(pUkf->predict.pFcnObserv[yIdx] != NULL)
            {
                //Propagate each sigma-point through observation Y_m[yL][2L+1] = h(X_m, u)             
                pUkf->predict.pFcnObserv[yIdx](&pUkf->input.u, &pUkf->predict.X_m, &pUkf->predict.Y_m,sigmaIdx);
            }
            else
            {
                //assign 0 if observation function is not specified
                pY_m[sigmaLen*sigmaIdx + yIdx]=0;
            }
            //Calculate mean of predicted output 
            //y_m[L] = sum(Wm(i)*Y_m(i))
            py_m[yIdx] += Wm[sigmaIdx] * pY_m[sigmaLen*sigmaIdx + yIdx];
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
                pPyy[yLen*yIdx + yTrIdx] += Wc[sigmaIdx]*term1*term2;            
            }           
        }
        
        //Calculate cross-covariance of state and output
        for(xIdx=0;xIdx<xLen;xIdx++)
        {
            for(yTrIdx=0;yTrIdx<yLen;yTrIdx++)
            {
                double term1 = (X_m[sigmaLen*xIdx + sigmaIdx] - px_m[xIdx]);
                double term2 = (pY_m[sigmaLen*yTrIdx + sigmaIdx] - py_m[yIdx]);

                
                Pxy[yLen*xIdx + yTrIdx] += Wc[sigmaIdx]*term1*term2;
            }
        }
    }
}

void ukf_meas_update(tUKF * const pUkf)
{
    tUKFupdate * pUpdate = (tUKFupdate*)&pUkf->update;


    //3.Calculate Kalman gain: 
    (void)mtx_identity_f64(&pUpdate->Ixx);
    
    //inv(Pyy)
    (void)mtx_inv_f64(&pUpdate->Pyy,&pUpdate->Ixx);

    //K = Pxy * inv(Pyy)
    (void)mtx_mul_f64(&pUpdate->Pxy,&pUpdate->Ixx, &pUpdate->K);
    
    //K'
    (void)mtx_transp_dest_f64(&pUpdate->K, &pUpdate->Kt);

    //4.Update state estimate

    // y = y - y_m
    (void)mtx_subtract_f64(&pUkf->input.y,&pUkf->predict.y_m);

    // K*(y - y_m) temporal stored in y
    (void)mtx_mul_f64(&pUpdate->K, &pUkf->input.y, &pUkf->input.y);

    // x_m + K*(y - y_m)
    (void)mtx_add_f64(&pUkf->predict.x_m, &pUkf->input.y);


    //5.Update error covariance
    //use Pxy for temporal resul trom multiplication
    //Pxy = K*Pyy
    (void)mtx_mul_f64(&pUpdate->K,&pUpdate->Pyy,&pUpdate->Pxy);

    //Ixx = K*Pyy*K'
    (void)mtx_mul_f64(&pUpdate->Pxy, &pUpdate->Kt, &pUpdate->Ixx);

    (void)mtx_subtract_f64(&pUkf->predict.X_m,&pUpdate->Ixx);

}




