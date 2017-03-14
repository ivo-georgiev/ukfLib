//This file depend of phisical system model

#include "ukfCfg.h"
#include "ukfLib.h"
#include "math.h"

tPredictFcn PredictFcn[stateVectorLen] = {&Fx0,&Fx1,&Fx2,&Fx3};
tObservFcn  ObservFcn[measVectorLen] = {&Hy1,&Hy2};

//-----------------------
//UKF Processing matrix
//-----------------------
double Wm_sigma_weight_1x9[9] = {3,3,3,3,3,3,3,3,3};
double Wc_sigma_weight_1x9[9] = {0,0,0,0,0,0,0,0,0};

//System input current u(k)
double u_curr_system_input_4x1[4] = {0,0,0,0}; 

//System input previous u(k)
double u_prev_system_input_4x1[4] = {0,0,0,0};

//System output measurement y(k)
double y_curr_system_meas_2x1[2] = {0,0};

//System output predicted y_m(k|k-1)
double y_mean_system_predict_2x1[2] = {0,0};

//System states: x(k), x(k-1), x(k|k-1) common array for all 
double x_system_states_4x1[4] = {0,0,50,50};

double x_system_states_correction_4x1[4] = {0,0,0,0};

//Sigma points X(k), X(k|k-1):
double X_sigma_points_4x9[4][9]=
{/*  s1  s2  s3  s4  s5  s6  s7  s8  s9        */
    {1,  2,  3,  4,  5,  6,  7,  8,  9}, /* x1 */
    {10,  11,  12,  13,  14,  15,  16,  17,  18}, /* x2 */
    {19,  20,  21,  22,  23,  24,  25,  26,  27}, /* x3 */
    {28,  29,  30,  31,  32,  33,  34,  35,  36}, /* x4 */
};

//Sigma points Y(k|k-1) = y_m
double Y_sigma_points_2x9[2][9]=
{/*  s1  s2  s3  s4  s5  s6  s7  s8  s9        */
    {0,  0,  0,  0,  0,  0,  0,  0,  0}, /* y1 */
    {0,  0,  0,  0,  0,  0,  0,  0,  0}, /* y2 */
};

//State covariance  P(k|k-1) = P_m, P(k)= P  
double Px_state_cov_4x4[4][4]=
{/*  x1, x2, x3, x4        */
    {0,  0,  0,  0}, /* x1 */
    {0,  0,  0,  0}, /* x2 */ 
    {0,  0,  0,  0}, /* x3 */  
    {0,  0,  0,  0}, /* x4 */
};

//State covariance initial values
double P0_state_cov_4x4[4][4]=
{/*  x1, x2, x3, x4        */
    {1,  0,  0,  0}, /* x1 */
    {0,  1,  0,  0}, /* x2 */ 
    {0,  0,  1,  0}, /* x3 */  
    {0,  0,  0,  1}, /* x4 */
};

//Process noise covariance Q : initial noise assumptions
double Qxx_process_noise_cov_4x4[4][4]=
{/*  x1, x2, x3, x4        */
    {0,  0,  0,  0}, /* x1 */
    {0,  0,  0,  0}, /* x2 */ 
    {0,  0,  4,  0}, /* x3 */  
    {0,  0,  0,  4}, /* x4 */
};

//Output noise covariance: initial noise assumptions
double Ryy_out_cov_noise_2x2[2][2]=
{/*  y1, y2         */
    {1,  0},  /* y1 */
    {0,  1},  /* y2 */
};

//Output covariance Pyy = R (initial assumption)
double Pyy_out_cov_2x2[2][2]=
{/*  y1, y2         */
    {0,  0},  /* y1 */
    {0,  0},  /* y2 */
};

double Pyy_out_cov_copy_2x2[2][2]=
{/*  y1, y2         */
    {0,  0},  /* y1 */
    {0,  0},  /* y2 */
};


//cross-covariance of state and output
double Pxy_state_out_cov_4x2[4][2]=
{/*  y1, y2         */
    {0,  0},  /* x1 */
    {0,  0},  /* x2 */
    {0,  0},  /* x3 */
    {0,  0},  /* x4 */
};

//Kalman gain matrix
double K_kalman_gain_4x2[4][2]=
{  
    {0, 0},
    {0, 0},
    {0, 0},
    {0, 0},
};

//Kalman gain transponce matrix
double K_kalman_transp_gain_2x4[2][4]=
{  
    {0, 0, 0, 0},
    {0, 0, 0, 0}
};

double Pxx_covariance_correction_4x4[4][4]=
{/*  x1, x2, x3, x4        */
    {0,  0,  0,  0}, /* x1 */
    {0,  0,  0,  0}, /* x2 */ 
    {0,  0,  0,  0}, /* x3 */  
    {0,  0,  0,  0}, /* x4 */
};

double temporal_2x2[2][2]=
{
    {0,  0},
    {0,  0},
};

//State transition matrix
static const double A[stateVectorLen][stateVectorLen] =
{
    {1,   0, dT0,   0},
    {0,   1,   0, dT0},
    {0,   0,   1,   0},
    {0,   0,   0,   1}
};



void Fx0(tMatrix * pu_p, tMatrix * pX_p, tMatrix * pX_m,int sigmaIdx)
{
    int nCol = pX_m->ncol; //pX_m->ncol == pX_p->ncol == 9

    //A[0][:]* X[:][0] - ToDo write function that mutiply specific rowXcol
    pX_m->val[nCol*0 + sigmaIdx] = (A[0][0] * pX_p->val[nCol*0 + sigmaIdx]) + (A[0][2] * pX_p->val[nCol*2 + sigmaIdx]);

}

void Fx1(tMatrix * pu_p, tMatrix * pX_p, tMatrix * pX_m,int sigmaIdx)
{
    int nCol = pX_m->ncol; //pX_m->ncol == pX_p->ncol == 9
   
    //A[1][:]* X[:][1]
    pX_m->val[nCol*1 + sigmaIdx] = (A[1][1] * pX_p->val[nCol*1 + sigmaIdx]) + (A[1][3] * pX_p->val[nCol*3 + sigmaIdx]);


}

void Fx2(tMatrix * pu_p, tMatrix * pX_p, tMatrix * pX_m,int sigmaIdx)
{
    int nCol = pX_m->ncol; //pX_m->ncol == pX_p->ncol == 9
    
    //fx2() = A[2][:]* X[:][2]
    pX_m->val[nCol*2 + sigmaIdx] = (A[2][2] * pX_p->val[nCol*2 + sigmaIdx]);
}

void Fx3(tMatrix * pu_p, tMatrix * pX_p, tMatrix * pX_m,int sigmaIdx)
{
    int nCol = pX_m->ncol; //pX_m->ncol == pX_p->ncol == 9
    
    //A[3][:]* X[:][3]
    pX_m->val[nCol*3 + sigmaIdx] = (A[3][3] * pX_p->val[nCol*3 + sigmaIdx]);
}

void Hy1(tMatrix * pu, tMatrix * pX_m, tMatrix * pY_m,int sigmaIdx)
{
    static const N1 = 20;
    static const E1 = 0;
    double term1;
    double term2;
    const int nCol = pY_m->ncol;

    term1 = pX_m->val[nCol*0 + sigmaIdx] - N1;
    term1 *= term1;

    term2 = pX_m->val[nCol*1 + sigmaIdx] - E1;
    term2 *= term2;

    pY_m->val[/*nCol*0 +*/ sigmaIdx] = sqrt(term1+term2);
    
}

void Hy2(tMatrix * pu, tMatrix * pX_m, tMatrix * pY_m,int sigmaIdx)
{
    static const N2 = 0;
    static const E2 = 20;
    double term1;
    double term2;
    const int nCol = pY_m->ncol;
    
    term1 = pX_m->val[nCol*0 + sigmaIdx] - N2;
    term1 *= term1;
    
    term2 = pX_m->val[nCol*1 + sigmaIdx] - E2;
    term2 *= term2;
    
    pY_m->val[nCol*1 + sigmaIdx] = sqrt(term1+term2);
    
}
