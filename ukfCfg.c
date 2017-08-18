//This file depend of phisical system model

#include "ukfCfg.h"
#include "ukfLib.h"
#include "math.h"

tPredictFcn PredictFcn[stateVectorLen] = {&Fx0,&Fx1,&Fx2,&Fx3};
tObservFcn  ObservFcn[measVectorLen] = {&Hy1,&Hy2};

//-----------------------
//UKF Processing matrix
//-----------------------
float64 Wm_sigma_weight_1x9[9] = {3,3,3,3,3,3,3,3,3};
float64 Wc_sigma_weight_1x9[9] = {0,0,0,0,0,0,0,0,0};

//System input current u(k)
float64 u_curr_system_input_4x1[4] = {0,0,0,0}; 

//System input previous u(k-1)
float64 u_prev_system_input_4x1[4] = {0,0,0,0};

//System output measurement y(k)
float64 y_curr_system_meas_2x1[2] = {0,0};

//System output predicted y_m(k|k-1)
float64 y_mean_system_predict_2x1[2] = {0,0};

//System states: x(k), x(k-1), x(k|k-1) common array for all 
double x_system_states_4x1[4] = {0,0,50,50};

float64 x_system_states_ic_4x1[4] = {0,0,50,50};

float64 x_system_states_correction_4x1[4] = {0,0,0,0};

//Sigma points X(k), X(k|k-1):
float64 X_sigma_points_4x9[4][9]=
{/*  s1  s2  s3  s4  s5  s6  s7  s8  s9        */
    {0,  0,  0,  0,  0,  0,  0,  0,  0}, /* x1 */
    {0,  0,  0,  0,  0,  0,  0,  0,  0}, /* x2 */
    {0,  0,  0,  0,  0,  0,  0,  0,  0}, /* x3 */
    {0,  0,  0,  0,  0,  0,  0,  0,  0}, /* x4 */
};

//Sigma points Y(k|k-1) = y_m
float64 Y_sigma_points_2x9[2][9]=
{/*  s1  s2  s3  s4  s5  s6  s7  s8  s9        */
    {0,  0,  0,  0,  0,  0,  0,  0,  0}, /* y1 */
    {0,  0,  0,  0,  0,  0,  0,  0,  0}, /* y2 */
};

//State covariance  P(k|k-1) = P_m, P(k)= P  
float64 Px_state_cov_4x4[4][4]=
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
float64 P0_state_cov_4x4[4][4]=
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
float64 Qxx_process_noise_cov_4x4[4][4]=
{/*  x1, x2, x3, x4        */
    {0,  0,  0,  0}, /* x1 */
    {0,  0,  0,  0}, /* x2 */ 
    {0,  0,  4,  0}, /* x3 */  
    {0,  0,  0,  4}, /* x4 */
};

//Output noise covariance: initial noise assumptions
float64 Ryy_out_cov_noise_2x2[2][2]=
{/*  y1, y2         */
    {1,  0},  /* y1 */
    {0,  1},  /* y2 */
};

//Output covariance Pyy = R (initial assumption)
float64 Pyy_out_cov_2x2[2][2]=
{/*  y1, y2         */
    {0,  0},  /* y1 */
    {0,  0},  /* y2 */
};

float64 Pyy_out_cov_copy_2x2[2][2]=
{/*  y1, y2         */
    {0,  0},  /* y1 */
    {0,  0},  /* y2 */
};


//cross-covariance of state and output
float64 Pxy_state_out_cov_4x2[4][2]=
{/*  y1, y2         */
    {0,  0},  /* x1 */
    {0,  0},  /* x2 */
    {0,  0},  /* x3 */
    {0,  0},  /* x4 */
};

//Kalman gain matrix
float64 K_kalman_gain_4x2[4][2]=
{  
    {0, 0},
    {0, 0},
    {0, 0},
    {0, 0},
};

//Kalman gain transponce matrix
float64 K_kalman_transp_gain_2x4[2][4]=
{  
    {0, 0, 0, 0},
    {0, 0, 0, 0}
};

float64 Pxx_covariance_correction_4x4[4][4]=
{/*  x1, x2, x3, x4        */
    {0,  0,  0,  0}, /* x1 */
    {0,  0,  0,  0}, /* x2 */ 
    {0,  0,  0,  0}, /* x3 */  
    {0,  0,  0,  0}, /* x4 */
};

float64 temporal_2x2[2][2]=
{
    {0,  0},
    {0,  0},
};

//State transition matrix
static const float64 A[stateVectorLen][stateVectorLen] =
{
    {1,   0, dT0,   0},
    {0,   1,   0, dT0},
    {0,   0,   1,   0},
    {0,   0,   0,   1}
};



void Fx0(tMatrix * pu_p, tMatrix * pX_p, tMatrix * pX_m,uint8 sigmaIdx)
{
    const uint8 nCol = pX_m->ncol; //pX_m->ncol == pX_p->ncol == 9

    //A[0][:]* X[:][0] - ToDo write function that mutiply specific rowXcol
    pX_m->val[nCol*0 + sigmaIdx] = (A[0][0] * pX_p->val[nCol*0 + sigmaIdx]) + (A[0][2] * pX_p->val[nCol*2 + sigmaIdx]);

}

void Fx1(tMatrix * pu_p, tMatrix * pX_p, tMatrix * pX_m,uint8 sigmaIdx)
{
    const uint8 nCol = pX_m->ncol; //pX_m->ncol == pX_p->ncol == 9
   
    //A[1][:]* X[:][1]
    pX_m->val[nCol*1 + sigmaIdx] = (A[1][1] * pX_p->val[nCol*1 + sigmaIdx]) + (A[1][3] * pX_p->val[nCol*3 + sigmaIdx]);


}

void Fx2(tMatrix * pu_p, tMatrix * pX_p, tMatrix * pX_m,uint8 sigmaIdx)
{
    const uint8 nCol = pX_m->ncol; //pX_m->ncol == pX_p->ncol == 9
    
    //fx2() = A[2][:]* X[:][2]
    pX_m->val[nCol*2 + sigmaIdx] = (A[2][2] * pX_p->val[nCol*2 + sigmaIdx]);
}

void Fx3(tMatrix * pu_p, tMatrix * pX_p, tMatrix * pX_m,uint8 sigmaIdx)
{
    const uint8 nCol = pX_m->ncol; //pX_m->ncol == pX_p->ncol == 9
    
    //A[3][:]* X[:][3]
    pX_m->val[nCol*3 + sigmaIdx] = (A[3][3] * pX_p->val[nCol*3 + sigmaIdx]);
}

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

    pY_m->val[/*nCol*0 +*/ sigmaIdx] = sqrt(term1+term2);
    
}

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
    
}
