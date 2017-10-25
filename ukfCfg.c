 /******************************************************************************************************************************************************************************************************\
 *** 
 *** Description       : IMPLEMENTATION OF THE ADDITIVE NOISE UKF: This example problem is "Computer Exercise 13.21" from (Simon, 2006).
 *** Codefile          : ukfCfg.c
 *** Documentation     : https://web.statler.wvu.edu/%7Eirl/IRL_WVU_Online_UKF_Implementation_V1.0_06_28_2013.pdf
 ***
 *** State vector
 *** state 0 = x[0] = n(k)
 *** state 1 = x[1] = e(k)
 *** state 2 = x[2] = ndot(k)
 *** state 3 = x[2] = edot(k)
 *** 
 *** 
 *** //State transition matrix
 *** {1,   0, dT0,   0},
 *** {0,   1,   0, dT0},
 *** {0,   0,   1,   0},
 *** {0,   0,   0,   1}
 *** 
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
#include "ukfCfg.h"

static void Fx1(tMatrix * pu_p, tMatrix * pX_p, tMatrix * pX_m,uint8 sigmaIdx, float64 dT);
static void Fx2(tMatrix * pu_p, tMatrix * pX_p, tMatrix * pX_m,uint8 sigmaIdx, float64 dT);
static void Fx3(tMatrix * pu_p, tMatrix * pX_p, tMatrix * pX_m,uint8 sigmaIdx, float64 dT);
static void Fx4(tMatrix * pu_p, tMatrix * pX_p, tMatrix * pX_m,uint8 sigmaIdx, float64 dT);

static void Hy1(tMatrix * pu, tMatrix * pX_m, tMatrix * pY_m,uint8 sigmaIdx);
static void Hy2(tMatrix * pu, tMatrix * pX_m, tMatrix * pY_m,uint8 sigmaIdx);


static tPredictFcn PredictFcn[4] = {&Fx1,&Fx2,&Fx3,&Fx4};
static tObservFcn  ObservFcn[2] = {&Hy1,&Hy2};

//-----------------------
//UKF Processing matrix
//-----------------------
static float64 Sc_vector[1][3] = {{1,2,0}};
static float64 Wm_weight_vector[1][9] = {{3,3,3,3,3,3,3,3,3}};
static float64 Wc_weight_vector[1][9] = {{0,0,0,0,0,0,0,0,0}};
static float64 u_system_input[4][1] = {{0},{0},{0},{0}}; 
static float64 u_prev_system_input[4][1] = {{0},{0},{0},{0}};
static float64 y_meas[2][1] = {{0},{0}};
static float64 y_predicted_mean[2][1] = {{0},{0}};
static float64 x_system_states[4][1] = {{0},{0},{50},{50}};
static float64 x_system_states_ic[4][1] = {{0},{0},{50},{50}};
static float64 x_system_states_correction[4][1] = {{0},{0},{0},{0}};
static float64 X_sigma_points[4][9]=
{/*  s1  s2  s3  s4  s5  s6  s7  s8  s9        */
    {0,  0,  0,  0,  0,  0,  0,  0,  0}, /* x1 */
    {0,  0,  0,  0,  0,  0,  0,  0,  0}, /* x2 */
    {0,  0,  0,  0,  0,  0,  0,  0,  0}, /* x3 */
    {0,  0,  0,  0,  0,  0,  0,  0,  0}, /* x4 */
};

//Sigma points Y(k|k-1) = y_m
static float64 Y_sigma_points[2][9]=
{/*  s1  s2  s3  s4  s5  s6  s7  s8  s9        */
    {0,  0,  0,  0,  0,  0,  0,  0,  0}, /* y1 */
    {0,  0,  0,  0,  0,  0,  0,  0,  0}, /* y2 */
};

//State covariance  P(k|k-1) = P_m, P(k)= P  
static float64 Pxx_error_covariance[4][4]=
{/*  x1, x2, x3, x4        */
    {0,  0,  0,  0}, /* x1 */
    {0,  0,  0,  0}, /* x2 */ 
    {0,  0,  0,  0}, /* x3 */  
    {0,  0,  0,  0}, /* x4 */
};

//State covariance initial values
/*Matthew B Rhudy : initial error covariance matrix should be defined based on your initialization error.
 I.e., if you think your initial state is not very close, the P0 value should be large,
 whereas if the initialization is very good (high confidence that your states are close to the correct values) you can assume a smaller P0 value.
  I would not recommend setting P0 to zero, as this assumes there is no initialization error and your initial states are perfect.  
  This is almost never the case.  Many times the P0 matrix is diagonal, with the diagonal components corresponding to the expected variance
  in the corresponding state, i.e. how much deviation you might expect in the initialization of that state.  If you have no idea where to start, 
 I recommend using an identity matrix rather than the zero matrix. */
static float64 Pxx0_init_error_covariance[4][4]=
{/*  x1, x2, x3, x4        */
    {1,  0,  0,  0}, /* x1 */
    {0,  1,  0,  0}, /* x2 */ 
    {0,  0,  1,  0}, /* x3 */  
    {0,  0,  0,  1}, /* x4 */
};

//Process noise covariance Q : initial noise assumptions
/* Matthew B Rhudy : Q matrix corresponds to the uncertainty that you expect in your state equations.
  This could include modeling errors or other uncertainties in the equations themselves.
  Some formulations consider input measurements in the state equations which introduces process noise.
  If you are very confident in your equations, you could set Q to zero. If you do that the filter will use
  the noise free model to predict the state vector and will ignore any measurement data since your model is assumed perfect. */
static float64 Qxx_process_noise_cov[4][4]=
{/*  x1, x2, x3, x4        */
    {0,  0,  0,  0}, /* x1 */
    {0,  0,  0,  0}, /* x2 */ 
    {0,  0,  4,  0}, /* x3 */  
    {0,  0,  0,  4}, /* x4 */
};

//Output noise covariance: initial noise assumptions
static float64 Ryy0_init_out_covariance[2][2]=
{/*  y1, y2         */
    {1,  0},  /* y1 */
    {0,  1},  /* y2 */
};

//Output covariance Pyy = R (initial assumption)
static float64 Pyy_out_covariance[2][2]=
{/*  y1, y2         */
    {0,  0},  /* y1 */
    {0,  0},  /* y2 */
};

static float64 Pyy_out_covariance_copy[2][2]=
{/*  y1, y2         */
    {0,  0},  /* y1 */
    {0,  0},  /* y2 */
};


//cross-covariance of state and output
static float64 Pxy_cross_covariance[4][2]=
{/*  y1, y2         */
    {0,  0},  /* x1 */
    {0,  0},  /* x2 */
    {0,  0},  /* x3 */
    {0,  0},  /* x4 */
};

//Kalman gain matrix
static float64 K_kalman_gain[4][2]=
{  
    {0, 0},
    {0, 0},
    {0, 0},
    {0, 0},
};

//Kalman gain transponce matrix
static float64 K_kalman_gain_transp[2][4]=
{  
    {0, 0, 0, 0},
    {0, 0, 0, 0}
};

static float64 _Pxx_covariance_correction_4x4[4][4]=
{/*  x1, x2, x3, x4        */
    {0,  0,  0,  0}, /* x1 */
    {0,  0,  0,  0}, /* x2 */ 
    {0,  0,  0,  0}, /* x3 */  
    {0,  0,  0,  0}, /* x4 */
};

static float64 I_identity_matrix[2][2]=
{
    {0,  0},
    {0,  0},
};

tUkfMatrix UkfMatrixCfg0 = 
{
    {NROWS(Sc_vector),NCOL(Sc_vector),&Sc_vector[0][0]},
    {NROWS(Wm_weight_vector),NCOL(Wm_weight_vector),&Wm_weight_vector[0][0]},
    {NROWS(Wc_weight_vector),NCOL(Wc_weight_vector),&Wc_weight_vector[0][0]},
    {NROWS(x_system_states),NCOL(x_system_states),&x_system_states[0][0]},
    {NROWS(x_system_states_ic),NCOL(x_system_states_ic),&x_system_states_ic[0][0]},
    {NROWS(x_system_states_correction),NCOL(x_system_states_correction),&x_system_states_correction[0][0]},
    {NROWS(u_system_input),NCOL(u_system_input),&u_system_input[0][0]},
    {NROWS(u_prev_system_input),NCOL(u_prev_system_input),&u_prev_system_input[0][0]},
    {NROWS(X_sigma_points),NCOL(X_sigma_points),&X_sigma_points[0][0]},
    {NROWS(Y_sigma_points),NCOL(Y_sigma_points),&Y_sigma_points[0][0]},
    {NROWS(y_predicted_mean),NCOL(y_predicted_mean),&y_predicted_mean[0][0]},   
    {NROWS(y_meas),NCOL(y_meas),&y_meas[0][0]},
    {NROWS(Pyy_out_covariance),NCOL(Pyy_out_covariance),&Pyy_out_covariance[0][0]},
    {NROWS(Pyy_out_covariance_copy),NCOL(Pyy_out_covariance_copy),&Pyy_out_covariance_copy[0][0]},
    {NROWS(Ryy0_init_out_covariance),NCOL(Ryy0_init_out_covariance),&Ryy0_init_out_covariance[0][0]},
    {NROWS(Pxy_cross_covariance),NCOL(Pxy_cross_covariance),&Pxy_cross_covariance[0][0]},
    {NROWS(Pxx_error_covariance),NCOL(Pxx_error_covariance),&Pxx_error_covariance[0][0]},
    {NROWS(Pxx0_init_error_covariance),NCOL(Pxx0_init_error_covariance),&Pxx0_init_error_covariance[0][0]},
    {NROWS(Qxx_process_noise_cov),NCOL(Qxx_process_noise_cov),&Qxx_process_noise_cov[0][0]},
    {NROWS(K_kalman_gain),NCOL(K_kalman_gain),&K_kalman_gain[0][0]},
    {NROWS(K_kalman_gain_transp),NCOL(K_kalman_gain_transp),&K_kalman_gain_transp[0][0]},
    {NROWS(I_identity_matrix),NCOL(I_identity_matrix),&I_identity_matrix[0][0]},  
    {NROWS(_Pxx_covariance_correction_4x4),NCOL(_Pxx_covariance_correction_4x4),&_Pxx_covariance_correction_4x4[0][0]},   
    &PredictFcn[0],
    &ObservFcn[0],
    0.1
};
/******************************************************************************************************************************************************************************************************\
 ***  FUNCTION:
 ***      void Fx1(tMatrix * pu_p, tMatrix * pX_p, tMatrix * pX_m,uint8 sigmaIdx)
 *** 
 ***  DESCRIPTION:
 ***       Calculate predicted state 0 for each sigma point. Note  that  this  problem  has  a  linear  prediction stage 
 ***       X_m[0][sigmaIdx] = f(X_p, u_p) = n(k)  = n(k-1) + dT * ndot(k-1)    
 ***            
 ***  PARAMETERS:
 ***      Type               Name              Range              Description
 ***      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 ***      tMatrix *          pu_p                                 Pointer to the input array/vector at (k-1) moment
 ***      tMatrix *          pX_p                                 Pointer to the sigma points array at (k-1) moment 
 ***      tMatrix *          pX_m                                 Pointer to the propagetad sigma points array at (k|k-1) moment (i.e prediction in moment k based on states in (k-1))
 ***      uint8              sigmaIdx                             Sigma point index.
 ***  RETURNS:
 ***           
 ***  SETTINGS:
 ***
\******************************************************************************************************************************************************************************************************/
void Fx1(tMatrix * pu_p, tMatrix * pX_p, tMatrix * pX_m,uint8 sigmaIdx, float64 dT)
{
    const uint8 nCol = pX_m->ncol; //pX_m->ncol == pX_p->ncol == 9

    pX_m->val[nCol*0 + sigmaIdx] = pX_p->val[nCol*0 + sigmaIdx] + dT * pX_p->val[nCol*2 + sigmaIdx];

	pu_p = pu_p; // todo: fix up
}
/******************************************************************************************************************************************************************************************************\
 ***  FUNCTION:
 ***      void Fx2(tMatrix * pu_p, tMatrix * pX_p, tMatrix * pX_m,uint8 sigmaIdx)
 *** 
 ***  DESCRIPTION:
 ***       Calculate predicted state 1 for each sigma point. Note  that  this  problem  has  a  linear  prediction stage 
 ***       X_m[1][sigmaIdx] = f(X_p, u_p) = e(k)  = e(k-1) + dT * edot(k-1)    
 ***            
 ***  PARAMETERS:
 ***      Type               Name              Range              Description
 ***      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 ***      tMatrix *          pu_p                                 Pointer to the input array/vector at (k-1) moment
 ***      tMatrix *          pX_p                                 Pointer to the sigma points array at (k-1) moment 
 ***      tMatrix *          pX_m                                 Pointer to the propagetad sigma points array at (k|k-1) moment (i.e prediction in moment k based on states in (k-1))
 ***      uint8              sigmaIdx                             Sigma point index.
 ***  RETURNS:
 ***           
 ***  SETTINGS:
 ***
\******************************************************************************************************************************************************************************************************/
void Fx2(tMatrix * pu_p, tMatrix * pX_p, tMatrix * pX_m,uint8 sigmaIdx, float64 dT)
{
    const uint8 nCol = pX_m->ncol; //pX_m->ncol == pX_p->ncol == 9
   
    pX_m->val[nCol*1 + sigmaIdx] = pX_p->val[nCol*1 + sigmaIdx] + dT * pX_p->val[nCol*3 + sigmaIdx];

	pu_p = pu_p; // todo: fix up
}
/******************************************************************************************************************************************************************************************************\
 ***  FUNCTION:
 ***      void Fx3(tMatrix * pu_p, tMatrix * pX_p, tMatrix * pX_m,uint8 sigmaIdx)
 *** 
 ***  DESCRIPTION:
 ***       Calculate predicted state 2 for each sigma point. Note  that  this  problem  has  a  linear  prediction stage 
 ***       X_m[2][sigmaIdx] = f(X_p, u_p) = ndot(k)  = ndot(k-1)    
 ***            
 ***  PARAMETERS:
 ***      Type               Name              Range              Description
 ***      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 ***      tMatrix *          pu_p                                 Pointer to the input array/vector at (k-1) moment
 ***      tMatrix *          pX_p                                 Pointer to the sigma points array at (k-1) moment 
 ***      tMatrix *          pX_m                                 Pointer to the propagetad sigma points array at (k|k-1) moment (i.e prediction in moment k based on states in (k-1))
 ***      uint8              sigmaIdx                             Sigma point index.
 ***  RETURNS:
 ***           
 ***  SETTINGS:
 ***
\******************************************************************************************************************************************************************************************************/
void Fx3(tMatrix * pu_p, tMatrix * pX_p, tMatrix * pX_m,uint8 sigmaIdx, float64 dT)
{
    const uint8 nCol = pX_m->ncol; //pX_m->ncol == pX_p->ncol == 9
    
    pX_m->val[nCol*2 + sigmaIdx] = pX_p->val[nCol*2 + sigmaIdx];
	pu_p = pu_p; // todo: fix up
}
/******************************************************************************************************************************************************************************************************\
 ***  FUNCTION:
 ***      void Fx4(tMatrix * pu_p, tMatrix * pX_p, tMatrix * pX_m,uint8 sigmaIdx)
 *** 
 ***  DESCRIPTION:
 ***       Calculate predicted state 3 for each sigma point. Note  that  this  problem  has  a  linear  prediction stage 
 ***       X_m[3][sigmaIdx] = f(X_p, u_p) = edot(k)  = edot(k-1)    
 ***            
 ***  PARAMETERS:
 ***      Type               Name              Range              Description
 ***      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 ***      tMatrix *          pu_p                                 Pointer to the input array/vector at (k-1) moment
 ***      tMatrix *          pX_p                                 Pointer to the sigma points array at (k-1) moment 
 ***      tMatrix *          pX_m                                 Pointer to the propagetad sigma points array at (k|k-1) moment (i.e prediction in moment k based on states in (k-1))
 ***      uint8              sigmaIdx                             Sigma point index.
 ***  RETURNS:
 ***           
 ***  SETTINGS:
 ***
\******************************************************************************************************************************************************************************************************/
void Fx4(tMatrix * pu_p, tMatrix * pX_p, tMatrix * pX_m,uint8 sigmaIdx, float64 dT)
{
    const uint8 nCol = pX_m->ncol; //pX_m->ncol == pX_p->ncol == 9
    
    pX_m->val[nCol*3 + sigmaIdx] = pX_p->val[nCol*3 + sigmaIdx];
	pu_p = pu_p; // todo: fix up
}
/******************************************************************************************************************************************************************************************************\
 ***  FUNCTION:
 ***      void Hy1(tMatrix * pu, tMatrix * pX_m, tMatrix * pY_m,uint8 sigmaIdx)
 *** 
 ***  DESCRIPTION:
 ***       Calculate predicted state 3 for each sigma point. This problem has a nonlinear observation 
 ***       Y_m[0][sigmaIdx] = h1(X_m, u) = y1(k)  = sqrt((n(k)-N1)^2+(e(k)-E1)^2)    
 ***            
 ***  PARAMETERS:
 ***      Type               Name              Range              Description
 ***      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 ***      tMatrix *          pu                                   Pointer to the input array/vector at (k) moment. Not used involved in system equation
 ***      tMatrix *          pY_m                                 Pointer to the predicted output at (k|k-1) moment (i.e prediction in moment k based on states in (k-1))
 ***      tMatrix *          pX_m                                 Pointer to the propagetad sigma points array at (k|k-1) moment (i.e prediction in moment k based on states in (k-1))
 ***      uint8              sigmaIdx                             Sigma point index.
 ***  RETURNS:
 ***           
 ***  SETTINGS:
 ***
\******************************************************************************************************************************************************************************************************/
void Hy1(tMatrix * pu, tMatrix * pX_m, tMatrix * pY_m,uint8 sigmaIdx)
{
    static const float64 N1 = 20;
    static const float64 E1 = 0;
    float64 term1;
    float64 term2;
    const uint8 nCol = pY_m->ncol;

    term1 = pX_m->val[nCol*0 + sigmaIdx] - N1;
    term1 *= term1;

    term2 = pX_m->val[nCol*1 + sigmaIdx] - E1;
    term2 *= term2;

    pY_m->val[sigmaIdx] = sqrt(term1+term2);
    
	pu = pu; // todo: fix up
}
/******************************************************************************************************************************************************************************************************\
 ***  FUNCTION:
 ***      void Hy2(tMatrix * pu, tMatrix * pX_m, tMatrix * pY_m,uint8 sigmaIdx)
 *** 
 ***  DESCRIPTION:
 ***       Calculate predicted state 3 for each sigma point. This problem has a nonlinear observation 
 ***       Y_m[1][sigmaIdx] = h2(X_m[0&1][sigmaIdx], u) = y2(k)  = sqrt( (n(k) - N2)^2 + (e(k) - E2)^2 )    
 ***            
 ***  PARAMETERS:
 ***      Type               Name              Range              Description
 ***      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 ***      tMatrix *          pu                                   Pointer to the input array/vector at (k) moment. Not used involved in system equation
 ***      tMatrix *          pY_m                                 Pointer to the predicted output at (k|k-1) moment (i.e prediction in moment k based on states in (k-1))
 ***      tMatrix *          pX_m                                 Pointer to the propagetad sigma points array at (k|k-1) moment (i.e prediction in moment k based on states in (k-1))
 ***      uint8              sigmaIdx                             Sigma point index.
 ***  RETURNS:
 ***           
 ***  SETTINGS:
 ***
\******************************************************************************************************************************************************************************************************/
void Hy2(tMatrix * pu, tMatrix * pX_m, tMatrix * pY_m,uint8 sigmaIdx)
{
    static const float64 N2 = 0;
    static const float64 E2 = 20;
    float64 term1;
    float64 term2;
    const uint8 nCol = pY_m->ncol;
    
    term1 = pX_m->val[nCol*0 + sigmaIdx] - N2;
    term1 *= term1;
    
    term2 = pX_m->val[nCol*1 + sigmaIdx] - E2;
    term2 *= term2;
    
    pY_m->val[nCol*1 + sigmaIdx] = sqrt(term1+term2);
    
	pu = pu; // todo: fix up
}
