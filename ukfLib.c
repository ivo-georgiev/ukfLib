/******************************************************************************************************************************************************************************************************\
 *** 
 *** Description       : IMPLEMENTATION OF THE ADDITIVE NOISE UKF
 *** Codefile          : ukfLib.c
 *** Documentation     : https://web.statler.wvu.edu/%7Eirl/IRL_WVU_Online_UKF_Implementation_V1.0_06_28_2013.pdf
 ***
 *** MIT License
 ***
 *** Copyright (c) 2017 ivo-georgiev
 ***  
 *** Permission is hereby granted, free of charge, to any person obtaining a copy
 *** of this software and associated documentation files (the "Software"), to deal
 *** in the Software without restriction, including without limitation the rights
 *** to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 *** copies of the Software, and to permit persons to whom the Software is
 *** furnished to do so, subject to the following conditions:
 ***    
 *** The above copyright notice and this permission notice shall be included in all
 *** copies or substantial portions of the Software.
 ***      
 *** THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 *** IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 *** FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 *** AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 *** LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 *** OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 *** SOFTWARE.
\******************************************************************************************************************************************************************************************************/

#include "ukfLib.h"

boolean ukf_dimension_check(tUKF * const pUkf);
boolean ukf_init(tUKF * const pUkf, tUkfMatrix * pUkfMatrix);
void ukf_step(tUKF * const pUkf);
void ukf_meas_update(tUKF * const pUkf);
void ukf_sigmapoint(tUKF * const pUkf);
void ukf_mean_pred_state(tUKF * const pUkf);
void ukf_mean_pred_output(tUKF * const pUkf);
void ukf_calc_covariances(tUKF * const pUkf);

/******************************************************************************************************************************************************************************************************\
 ***  FUNCTION:
 ***      boolean ukf_dimension_check(tUKF * const pUkf)
 *** 
 ***  DESCRIPTION:
 ***      Check if working matrix size defined in ukfCfg.c match to defined system expectation(verification of all matrix is vs number of states xLen and number of measurements yLen)     
 ***            
 ***  PARAMETERS:
 ***      Type               Name              Range              Description
 ***      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 ***      tUKF * const       pUkf                                 UKF - Working structure with reference to all in,out,states,par
 ***  RETURNS:
 ***      boolean
 ***      0 - OK , 
 ***      1 - NOK
 ***  SETTINGS:
 ***
\******************************************************************************************************************************************************************************************************/
boolean ukf_dimension_check(tUKF * const pUkf)
{
    const uint8  stateLen = pUkf->par.xLen;
    const uint8  measLen = pUkf->par.yLen;
    const uint8  sigmaLen = pUkf->par.sLen;
    boolean Result = 0;

    //check system input vector size: (xLen x 1)
    if((pUkf->input.u.nrow != stateLen || pUkf->input.u.ncol != 1)&&
       (pUkf->prev.u_p.nrow != stateLen || pUkf->prev.u_p.ncol != 1))
    {
        Result |= 1;
    } 
    //check measurement vector size: (yLen x 1)
    if(pUkf->input.y.nrow != pUkf->par.yLen || pUkf->input.y.ncol != 1)
    {
        Result |= 1;
    }  
    //check Wm,Wc sigma weight matrix size: (1 x sLen)
    if((pUkf->par.Wm.nrow != 1 || pUkf->par.Wm.ncol != sigmaLen) &&
       (pUkf->par.Wc.nrow != 1 || pUkf->par.Wc.ncol != sigmaLen))
    {
        Result |= 1;
    }  
    //check initial covariance matrix size: (xLen x xLen)
    if(pUkf->par.Pxx0.nrow != stateLen || pUkf->par.Pxx0.ncol != stateLen)
    {
        Result |= 1;
    }
    //check Process noise covariance Q: (xLen x xLen)
    if(pUkf->par.Qxx.nrow != stateLen || pUkf->par.Qxx.ncol != stateLen)
    {
        Result |= 1;
    }
    //check Output noise covariance matrix size: (yLen x yLen)
    if(pUkf->par.Ryy0.nrow != pUkf->par.yLen || pUkf->par.Ryy0.ncol != pUkf->par.yLen)
    {
        Result |= 1;
    }
    //check X sigma point matrix size: (xLen x 2*xLen+1)
    if(pUkf->predict.X_m.nrow != stateLen || pUkf->predict.X_m.ncol != pUkf->par.sLen)
    {
        Result |= 1;
    }
    //check Y sigma point matrix size: (yLen x 2*xLen+1) , Y(k|k-1) = y_m
    if(pUkf->predict.Y_m.nrow != pUkf->par.yLen || pUkf->predict.Y_m.ncol != pUkf->par.sLen)
    {
        Result |= 1;
    }
    //check state/error covariance matrix size: (xLen x xLen) , Pxx_p == P_m == Pxx
    if(pUkf->predict.P_m.nrow != stateLen || pUkf->predict.P_m.ncol != stateLen)
    {
        Result |= 1;
    }
    //check Output covariance and it's copy size 
    if((pUkf->update.Pyy.nrow != pUkf->par.yLen || pUkf->update.Pyy.ncol != pUkf->par.yLen) &&
       (pUkf->update.Pyy_cpy.nrow != pUkf->par.yLen || pUkf->update.Pyy_cpy.ncol != pUkf->par.yLen))
    {
        Result |= 1;
    }
    //check cross-covariance matrix of state and output size: (xLen x yLen)
    if(pUkf->update.Pxy.nrow != stateLen || pUkf->update.Pxy.ncol != pUkf->par.yLen)
    {
        Result |= 1;
    }
    //check Pxx covariance correction: (xLen x xLen)
    if(pUkf->update.Pxx_corr.nrow != stateLen || pUkf->update.Pxx_corr.ncol != stateLen)
    {
        Result |= 1;
    }
    //check Kalman gain matrix and it's transp: (xLen x yLen)
    if((pUkf->update.K.nrow != stateLen || pUkf->update.K.ncol != pUkf->par.yLen) &&
       (pUkf->update.Kt.nrow != pUkf->par.yLen || pUkf->update.Kt.ncol != stateLen))
    {
        Result |= 1;
    }

    return Result;
}
/******************************************************************************************************************************************************************************************************\
 ***  FUNCTION:
 ***      ukf_init
 *** 
 ***  DESCRIPTION:
 ***      Unscented Kalman filter initialization. Calculate some filter parameter and initialize matrix pointers with real working arrays. Check matrix/array size consistency  
 *** 
 ***  PARAMETERS:
 ***      Type               Name                     Lim       Description
 ***      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 ***      tUKF * const       pUkf                               UKF - Working structure with all in,out,par
 ***      uint8              xLen                     1-127     UKF - State vector length
 ***      uint8              yLen                               UKF - Measurements vector length
 ***      tUkfMatrix *       pUkfMatrix                         UKF - Structure with all filter matrix
 ***  RETURNS:
 ***      boolean
 ***
 ***  SETTINGS:
 ***
\******************************************************************************************************************************************************************************************************/
boolean ukf_init(tUKF * const pUkf, tUkfMatrix * pUkfMatrix)
{
    tUKFpar * const pPar = (tUKFpar *)&pUkf->par;
    tUKFprev * const pPrev = (tUKFprev *)&pUkf->prev;
    const uint8 WmLen = pUkfMatrix->Wm_weight_vector.ncol;
    const uint8 WcLen = pUkfMatrix->Wc_weight_vector.ncol;
    
    pPar->x0 = pUkfMatrix->x_system_states_ic;
    pPar->Ryy0 = pUkfMatrix->Ryy0_init_out_covariance;
    pPar->Pxx0 = pUkfMatrix->Pxx0_init_error_covariance;
    pPar->Qxx = pUkfMatrix->Qxx_process_noise_cov;
    pPar->Wm =  pUkfMatrix->Wm_weight_vector;
    pPar->Wc =  pUkfMatrix->Wc_weight_vector;
    pPar->alpha = pUkfMatrix->Sc_vector.val[alphaIdx];//scaling[alphaIdx];
    pPar->betha = pUkfMatrix->Sc_vector.val[bethaIdx];;//scaling[bethaIdx];
    pPar->kappa = pUkfMatrix->Sc_vector.val[kappaIdx];;//scaling[kappaIdx];
    pPar->xLen = pUkfMatrix->x_system_states.nrow; 
    pPar->yLen = pUkfMatrix->y_predicted_mean.nrow;
    pPar->sLen = 2*pPar->xLen+1;
    
    //#1.3'(begin) Calculate scaling parameter
    pPar->lambda = pPar->alpha * pPar->alpha;
    pPar->lambda *= (float64)(pPar->xLen + pPar->kappa);
    pPar->lambda -= (float64)pPar->xLen;
    //#1.3'(end) Calculate scaling parameter
    
    //#1.2'(begin) Calculate weight vectors
    if(WmLen == pPar->sLen && WcLen == WmLen)
    {
        uint8 col;
        const float64 Wm0 = pPar->lambda/(pPar->xLen + pPar->lambda);

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
    //#1.2'(end) Calculate weight vectors

    pUkf->input.u = pUkfMatrix->u_system_input;
    pUkf->input.y = pUkfMatrix->y_meas;

    pPrev->Pxx_p = pUkfMatrix->Pxx_error_covariance;
    pPrev->X_p = pUkfMatrix->X_sigma_points; //share same memory with X_m
    pPrev->u_p = pUkfMatrix->u_system_input;//u_prev_system_input;
    pPrev->x_p = pUkfMatrix->x_system_states;

    pUkf->predict.P_m = pUkfMatrix->Pxx_error_covariance;
    pUkf->predict.X_m = pUkfMatrix->X_sigma_points;
    pUkf->predict.x_m = pUkfMatrix->x_system_states;
    pUkf->predict.Y_m = pUkfMatrix->Y_sigma_points;
    pUkf->predict.y_m = pUkfMatrix->y_predicted_mean;
    pUkf->predict.pFcnPredict = pUkfMatrix->fcnPredict;
    pUkf->predict.pFcnObserv = pUkfMatrix->fcnObserve;

    pUkf->update.Iyy = pUkfMatrix->I_identity_matrix;
    pUkf->update.K = pUkfMatrix->K_kalman_gain;
    pUkf->update.Kt = pUkfMatrix->K_kalman_gain_transp;
    pUkf->update.Pxx = pUkfMatrix->Pxx_error_covariance;
    pUkf->update.Pxy = pUkfMatrix->Pxy_cross_covariance;
    pUkf->update.Pyy = pUkfMatrix->Pyy_out_covariance;
    pUkf->update.Pyy_cpy = pUkfMatrix->Pyy_out_covariance_copy;
    pUkf->update.x = pUkfMatrix->x_system_states; //&px = &px_m = &px_p
    pUkf->update.x_corr = pUkfMatrix->x_system_states_correction;
    pUkf->update.Pxx_corr = pUkfMatrix->Pxx_covariance_correction;

    mtx_cpy_f64(&pUkf->prev.Pxx_p, &pPar->Pxx0);//init also P_m, Pxx   
    //mtx_cpy_f64(&pUkf->prev.x_p, &pPar->x0);//init also x_m
    /*mtx_cpy_f64(&pUkf->update.Pyy, &pPar->Ryy0);
    mtx_zeros_f64(&pUkf->prev.X_p);//inti also X_m
    mtx_zeros_f64(&pUkf->prev.u_p); 
    mtx_zeros_f64(&pUkf->predict.y_m);
    mtx_zeros_f64(&pUkf->predict.Y_m);
    mtx_zeros_f64(&pUkf->update.Iyy);
    mtx_zeros_f64(&pUkf->update.K);
    mtx_zeros_f64(&pUkf->update.Kt);
    mtx_zeros_f64(&pUkf->update.Pxy);  
    mtx_zeros_f64(&pUkf->update.Pyy_cpy);
    mtx_zeros_f64(&pUkf->update.x_corr);
    mtx_zeros_f64(&pUkf->update.Pxx_corr);     
    mtx_zeros_f64(&pUkf->input.u);
    mtx_zeros_f64(&pUkf->input.y);*/

    return ukf_dimension_check(pUkf);  
}
/******************************************************************************************************************************************************************************************************\
 ***  FUNCTION:
 ***      ukf_step(tUKF * const pUkf)
 *** 
 ***  DESCRIPTION:
 ***      Unscented Kalman filter periodic task. All new inputs(measurements, system inputs) should be updated periodicaly before execution of ukf_step().
 ***      UKF processing is separated on two sub-steps
 ***      - Predict
 ***      - Measurement update
 ***  PARAMETERS:
 ***      Type               Name              Range              Description
 ***      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 ***      tUKF * const       pUkf                                 UKF - Working structure with reference to all in,out,states,par
 ***  RETURNS:
 ***      void
 ***
 ***  SETTINGS:
 ***
\******************************************************************************************************************************************************************************************************/
void ukf_step(tUKF * const pUkf)
{
    float64 * const pu_p = pUkf->prev.u_p.val;
    const float64 * const pu = pUkf->input.u.val;
    const uint8 uLen = pUkf->prev.u_p.nrow;//input length
    uint8 u8Idx;

    ukf_sigmapoint(pUkf);
    ukf_mean_pred_state(pUkf);
    ukf_mean_pred_output(pUkf);
    ukf_calc_covariances(pUkf);
    ukf_meas_update(pUkf);
    
    for(u8Idx=0;u8Idx<uLen;u8Idx++)
    {
        //store prev inputs required for next step calculation
        pu_p[u8Idx] = pu[u8Idx];
    }
}

/******************************************************************************************************************************************************************************************************\
 ***  FUNCTION:
 ***      void ukf_sigmapoint(tUKF * const pUkf)
 *** 
 ***  DESCRIPTION:
 ***      Step 1:  Generate the Sigma-Points
 ***      # 1.1 Calculate error covariance matrix square root : sqrt(Pxx_p) = chol(Pxx_p) 
 ***      # 1.2 Calculate the sigma-points                    : X_p[L][2L+1] == X(k-1) ,wher L is number of system states xLen      
 ***            
 ***  PARAMETERS:
 ***      Type               Name              Range              Description
 ***      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 ***      tUKF * const       pUkf                                 UKF - Working structure with reference to all in,out,states,par
 ***  RETURNS:
 ***      void
 ***
 ***  SETTINGS:
 ***
\******************************************************************************************************************************************************************************************************/
void ukf_sigmapoint(tUKF * const pUkf)
{
    float64 * const pPxx_p = pUkf->prev.Pxx_p.val;
    float64 * const pX_p = pUkf->prev.X_p.val;
    float64 * const px_p = pUkf->prev.x_p.val;
    float64 f64Scalar;
    const float64 lambda = pUkf->par.lambda;
    const uint8 sLen = pUkf->par.sLen;
    const uint8 xLen = pUkf->par.xLen;
    uint8 xIdx;
    uint8 sigmaIdx=0;
    mtxResultInfo mtxResult;
    
    //#1.1(begin) Calculate error covariance matrix square root
    mtxResult = mtx_chol_lower_f64(&pUkf->prev.Pxx_p);

    //#1.1(end) Calculate error covariance matrix square root

    if(MTX_OPERATION_OK == mtxResult)
    {       
        //#1.2(begin) Calculate the sigma-points
        for(xIdx=0;xIdx<xLen;xIdx++)
        {
            //first column of sigma point matrix is equal of previous state value
            pX_p[sLen*xIdx+sigmaIdx] = px_p[xIdx];
        }
               
        f64Scalar = sqrt(xLen + lambda);
        
        (void)mtx_mul_scalar_f64(&pUkf->prev.Pxx_p,f64Scalar);
        
        for(sigmaIdx=1;sigmaIdx < sLen;sigmaIdx++)
        {
            for(xIdx=0;xIdx<xLen;xIdx++)
            {
                if(sigmaIdx <= xLen)
                {
                    pX_p[sLen*xIdx+sigmaIdx] = px_p[xIdx] + pPxx_p[xLen*xIdx + (sigmaIdx-1)];
                }
                else
                {
                    pX_p[sLen*xIdx+sigmaIdx] = px_p[xIdx] - pPxx_p[xLen*xIdx + (sigmaIdx-xLen-1)];
                }            
            }
        }
        //#1.2(end) Calculate the sigma-points
    }
    else
    {
    }
}
/******************************************************************************************************************************************************************************************************\
 ***  FUNCTION:
 ***      void ukf_mean_pred_state(tUKF * const pUkf)
 *** 
 ***  DESCRIPTION:
 ***      Step 2: Prediction Transformation (APPENDIX A:IMPLEMENTATION OF THE ADDITIVE NOISE UKF)
 ***      # 2.1 Propagate each sigma-point through prediction  : X_m = f(X_p, u_p)
 ***      # 2.2 Calculate mean of predicted state              : x_m = sum(Wm(i)*X_m(i)) , i=0,..2L
 ***            
 ***  PARAMETERS:
 ***      Type               Name              Range              Description
 ***      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 ***      tUKF * const       pUkf                                 UKF - Working structure with reference to all in,out,states,par
 ***  RETURNS:
 ***      void
 ***  SETTINGS:
 ***
\******************************************************************************************************************************************************************************************************/ 
void ukf_mean_pred_state(tUKF * const pUkf)
{
    tUKFpar const * const pPar = (tUKFpar *)&pUkf->par;
    const uint8 xLen = pPar->xLen;
    const uint8 sigmaLen = pPar->sLen;
    float64 * const px_m = pUkf->predict.x_m.val;
    float64 const * const pX_m = pUkf->predict.X_m.val;
    float64 const * const pWm = pPar->Wm.val;
    uint8 sigmaIdx,xIdx;

    for(xIdx=0;xIdx<xLen;xIdx++)
    {
        px_m[xIdx] = 0;
        
        for(sigmaIdx=0;sigmaIdx<sigmaLen;sigmaIdx++)
        {                     
            if(pUkf->predict.pFcnPredict[xIdx] != NULL)
            {
                //#2.1 Propagate each sigma-point through prediction
                pUkf->predict.pFcnPredict[xIdx](&pUkf->prev.u_p, &pUkf->prev.X_p, &pUkf->predict.X_m,sigmaIdx);
            }
            //#2.2 Calculate mean of predicted state 
            px_m[xIdx] += pWm[sigmaIdx] * pX_m[sigmaLen*xIdx + sigmaIdx];
        }
    }
}
/******************************************************************************************************************************************************************************************************\
 ***  FUNCTION:
 ***      void ukf_mean_pred_output(tUKF * const pUkf)
 *** 
 ***  DESCRIPTION:
 ***      Step 3: Observation Transformation (APPENDIX A:IMPLEMENTATION OF THE ADDITIVE NOISE UKF)
 ***      # 2.3 Calculate covariance of predicted state        : P_m = Wc(sigmaIdx)*(X_m-x_m)*(X_m-x_m)'    P(k|k-1)
 ***      # 3.1 Propagate each sigma-point through observation : Y_m = h(X_m, u)
 ***      # 3.2 Calculate mean of predicted output             : y_m = sum(Wm(i)*Y_m(i))
 ***            
 ***  PARAMETERS:
 ***      Type               Name              Range              Description
 ***      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 ***      tUKF * const       pUkf                                 UKF - Working structure with reference to all in,out,states,par
 ***  RETURNS:
 ***      void
 ***  SETTINGS:
 ***
\******************************************************************************************************************************************************************************************************/
void ukf_mean_pred_output(tUKF * const pUkf)
{
    tUKFpar const * const pPar = (tUKFpar *)&pUkf->par;
    float64 const * const pWm = pPar->Wm.val;
    float64 const * const pWc = pPar->Wc.val;
    float64 const * const pX_m = pUkf->predict.X_m.val;
    float64 * const pY_m = pUkf->predict.Y_m.val;
    float64 * const pP_m = pUkf->predict.P_m.val;
    float64 * const px_m = pUkf->predict.x_m.val;
    float64 * py_m = pUkf->predict.y_m.val;
    const uint8 sigmaLen = pPar->sLen;
    const uint8 xLen = pPar->xLen;
    const uint8 yLen= pPar->yLen;
    uint8 sigmaIdx,xIdx,xTrIdx,yIdx;
 
    memset(py_m,0,sizeof(float64)*yLen); 
    
    //P(k|k-1) = Q(k-1)
    mtx_cpy_f64(&pUkf->predict.P_m, &pUkf->par.Qxx);
    
    for(sigmaIdx=0;sigmaIdx<sigmaLen;sigmaIdx++)
    {             
        for(xIdx=0;xIdx<xLen;xIdx++)
        {                                 
            for(xTrIdx=0;xTrIdx<xLen;xTrIdx++)
            { 
                float64 term1 = (pX_m[sigmaLen*xIdx + sigmaIdx] - px_m[xIdx]);
                float64 term2 = (pX_m[sigmaLen*xTrIdx + sigmaIdx] - px_m[xTrIdx]);
                
                //#2.3 Calculate covariance of predicted state
                //Perform multiplication with accumulation for each covariance matrix index
                pP_m[xLen*xIdx + xTrIdx] += pWc[sigmaIdx]*term1*term2;
            }           
        } 
        
        for(yIdx=0;yIdx<yLen;yIdx++)
        {
            if(pUkf->predict.pFcnObserv[yIdx] != NULL)
            {
                //#3.1 Propagate each sigma-point through observation             
                pUkf->predict.pFcnObserv[yIdx](&pUkf->input.u, &pUkf->predict.X_m, &pUkf->predict.Y_m,sigmaIdx);
            }
            else
            {
                //assign 0 if observation function is not specified
                pY_m[sigmaLen*sigmaIdx + yIdx]=0;
            }
            //#3.2 Calculate mean of predicted output 
            py_m[yIdx] += pWm[sigmaIdx] * pY_m[sigmaLen*yIdx + sigmaIdx];
        }
    }
}
/******************************************************************************************************************************************************************************************************\
 ***  FUNCTION:
 ***      void ukf_calc_covariances(tUKF * const pUkf)
 *** 
 ***  DESCRIPTION:
 ***      # 3.3 Calculate covariance of predicted output       : Pyy = Wc(sigmaIdx)*(Y_m-y_m)*(Y_m-y_m)'
 ***      # 3.4 Calculate cross-covariance of state and output : Pxy = Q + sum(Wc*()*()')
 ***            
 ***  PARAMETERS:
 ***      Type               Name              Range              Description
 ***      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 ***      tUKF * const       pUkf                                 UKF - Working structure with reference to all in,out,states,par
 ***  RETURNS:
 ***      void
 ***  SETTINGS:
 ***
\******************************************************************************************************************************************************************************************************/ 
void ukf_calc_covariances(tUKF * const pUkf)
{
    tUKFpar const * const pPar = (tUKFpar *)&pUkf->par;
    float64 const * const pWc = pPar->Wc.val;
    float64 const * const pX_m = pUkf->predict.X_m.val;
    float64 * const pY_m = pUkf->predict.Y_m.val;
    float64 * const pPyy = pUkf->update.Pyy.val;
    float64 * const pPxy = pUkf->update.Pxy.val;
    float64 * const px_m = pUkf->predict.x_m.val;
    float64 * py_m = pUkf->predict.y_m.val;
    const uint8 sigmaLen = pPar->sLen;
    const uint8 xLen = pPar->xLen;
    const uint8 yLen= pPar->yLen;
    uint8 sigmaIdx,xIdx,yIdx,yTrIdx;

    mtx_cpy_f64(&pUkf->update.Pyy, &pPar->Ryy0);//Pyy(k|k-1) = R(k)
    
    memset(pPxy,0,sizeof(float64)*xLen*yLen);
    
    for(sigmaIdx=0;sigmaIdx<sigmaLen;sigmaIdx++)
    {              
        for(yIdx=0;yIdx<yLen;yIdx++)
        {                      
            for(yTrIdx=0;yTrIdx<yLen;yTrIdx++)
            {
                //loop col of COV[L][:]                
                float64 term1 = (pY_m[sigmaLen*yIdx + sigmaIdx] - py_m[yIdx]);
                float64 term2 = (pY_m[sigmaLen*yTrIdx + sigmaIdx] - py_m[yTrIdx]);
                
                //#3.3 Calculate covariance of predicted output
                //Perform multiplication with accumulation for each covariance matrix index 
                pPyy[yLen*yIdx + yTrIdx] += pWc[sigmaIdx]*term1*term2;            
            }           
        }
        
        for(xIdx=0;xIdx<xLen;xIdx++)
        {
            for(yTrIdx=0;yTrIdx<yLen;yTrIdx++)
            {
                float64 term1 = (pX_m[sigmaLen*xIdx + sigmaIdx] - px_m[xIdx]);
                float64 term2 = (pY_m[sigmaLen*yTrIdx + sigmaIdx] - py_m[yTrIdx]);
                
                //#3.4 Calculate cross-covariance of state and output
                pPxy[yLen*xIdx + yTrIdx] += pWc[sigmaIdx]*term1*term2;
            }
        }
    }
}
/******************************************************************************************************************************************************************************************************\
 ***  FUNCTION:
 ***      void ukf_meas_update(tUKF * const pUkf)
 *** 
 ***  DESCRIPTION:
 ***      Step 4: Measurement Update (APPENDIX A:IMPLEMENTATION OF THE ADDITIVE NOISE UKF)
 ***      # 4.1 Calculate Kalman gain          : K = Pxy*inv(Pyy)
 ***      # 4.2 Update state estimate          : x = x_m + K(y - y_m)
 ***      # 4.3 Update error covariance        : Pxx = Pxx_m - K*Pyy*K'
 ***            
 ***  PARAMETERS:
 ***      Type               Name              Range              Description
 ***      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 ***      tUKF * const       pUkf                                 UKF - Working structure with reference to all in,out,states,par
 ***  RETURNS:
 ***      void
 ***  SETTINGS:
 ***
\******************************************************************************************************************************************************************************************************/
void ukf_meas_update(tUKF * const pUkf)
{
    tUKFupdate * const pUpdate = (tUKFupdate*)&pUkf->update;

    //#4.1(begin) Calculate Kalman gain:
    (void)mtx_identity_f64(&pUpdate->Iyy);
    
    (void)mtx_cpy_f64(&pUpdate->Pyy_cpy, &pUpdate->Pyy);

    //inv(Pyy_cpy)
    (void)mtx_inv_f64(&pUpdate->Pyy_cpy,&pUpdate->Iyy);

    //Kgain = Pxy * inv(Pyy)
    (void)mtx_mul_f64(&pUpdate->Pxy,&pUpdate->Iyy, &pUpdate->K);
    
    //K' - transponed gain needed for further calculationd
    (void)mtx_transp_dest_f64(&pUpdate->K, &pUpdate->Kt);
    //#4.1(end) Calculate Kalman gain:

    //#4.2(begin) Update state estimate
    // y = y - y_m
    (void)mtx_subtract_f64(&pUkf->input.y,&pUkf->predict.y_m);

    // K*(y - y_m) states correction 
    (void)mtx_mul_f64(&pUpdate->K, &pUkf->input.y, &pUkf->update.x_corr);

    // x = x_m + K*(y - y_m)
    (void)mtx_add_f64(&pUkf->predict.x_m, &pUkf->update.x_corr);
    //#4.2(end) Update state estimate

    //#4.3(begin).Update error covariance
    //use Pxy for temporal result from multiplication
    //Pxy = K*Pyy
    (void)mtx_mul_f64(&pUpdate->K,&pUpdate->Pyy,&pUpdate->Pxy);

    //Pxx_corr = K*Pyy*K'
    (void)mtx_mul_f64(&pUpdate->Pxy, &pUpdate->Kt, &pUpdate->Pxx_corr);

    (void)mtx_subtract_f64(&pUkf->predict.P_m,&pUpdate->Pxx_corr);
    //#4.3(end).Update error covariance
}




