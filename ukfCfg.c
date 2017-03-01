//This file depend of phisical system model

#include "ukfCfg.h"
#include "ukfLib.h"
#include "math.h"

void Fx0(tMatrix * pu_p, tMatrix * pX_p, tMatrix * pX_m,int sigmaIdx);
void Fx1(tMatrix * pu_p, tMatrix * pX_p, tMatrix * pX_m,int sigmaIdx);
void Fx2(tMatrix * pu_p, tMatrix * pX_p, tMatrix * pX_m,int sigmaIdx);
void Fx3(tMatrix * pu_p, tMatrix * pX_p, tMatrix * pX_m,int sigmaIdx);

void Hy1(tMatrix * pu, tMatrix * pX_m, tMatrix * pY_m,int sigmaIdx);
void Hy2(tMatrix * pu, tMatrix * pX_m, tMatrix * pY_m,int sigmaIdx);

static const tPredictFcn PredictFcn[stateVectorLen] = {&Fx0,&Fx1,&Fx2,&Fx3};
static const tObservFcn  ObservFcn[stateVectorLen] = {&Hy1,&Hy2};

//-----------------------
//UKF Processing matrix
//-----------------------
double sigma_weight_Wm_1x9[1][9] = {0,0,0,0,0,0,0,0,0};
double sigma_weight_Wc_1x9[1][9] = {0,0,0,0,0,0,0,0,0};

//System input current u(k)
double u_curr_system_input_4x1[4] = {0,0,0,0}; 

//System input previous u(k)
double u_prev_system_input_4x1[4] = {0,0,0,0};

//System output measurement y(k)
double y_curr_system_meas_2x1[4] = {0,0};

//System output predicted y_m(k|k-1)
double y_mean_system_predict_2x1[4] = {0,0};

//System states: x(k), x(k-1), x(k|k-1) common array for all 
double x_system_states_4x1[4] = {0,0,0,0};

//Sigma points X(k), X(k|k-1):
double X_sigma_points_4x9[4][9]=
{/*  s1  s2  s3  s4  s5  s6  s7  s8  s9        */
    {0,  0,  0,  0,  0,  0,  0,  0,  0}, /* x1 */
    {0,  0,  0,  0,  0,  0,  0,  0,  0}, /* x2 */
    {0,  0,  0,  0,  0,  0,  0,  0,  0}, /* x3 */
    {0,  0,  0,  0,  0,  0,  0,  0,  0}, /* x4 */
};

//Sigma points Y(k|k-1) = y_m
double Y_sigma_points_2x9[4][9]=
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
const double P0_state_cov_4x4[4][4]=
{/*  x1, x2, x3, x4        */
    {1,  0,  0,  0}, /* x1 */
    {0,  1,  0,  0}, /* x2 */ 
    {0,  0,  1,  0}, /* x3 */  
    {0,  0,  0,  1}, /* x4 */
};

//Process noise covariance Q : initial noise assumptions
const double Qx_process_noise_cov_4x4[4][4]=
{/*  x1, x2, x3, x4        */
    {0,  0,  0,  0}, /* x1 */
    {0,  0,  0,  0}, /* x2 */ 
    {0,  0,  4,  0}, /* x3 */  
    {0,  0,  0,  4}, /* x4 */
};

//Output noise covariance: initial noise assumptions
const double Ry_out_cov_noise_2x2[2][2]=
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

//temporal result ,matrix 4*4
double temp_4x4[4][4]=
{
    {0, 0, 0, 0},
    {0, 0, 0, 0},
    {0, 0, 0, 0},
    {0, 0, 0, 0},
};

//temporal result ,matrix 4*4
const double Identity_4x4[4][4]=
{
    {1,0,0,0},
    {0,1,0,0},
    {0,0,1,0},
    {0,0,0,1},
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
    pX_m->val[nCol*sigmaIdx + idxSt0] = (A[0][0] * pX_p->val[nCol*idxSt0 + sigmaIdx]) + (A[0][2] * pX_p->val[nCol*idxSt2 + sigmaIdx]);

}

void Fx1(tMatrix * pu_p, tMatrix * pX_p, tMatrix * pX_m,int sigmaIdx)
{
    int nCol = pX_m->ncol; //pX_m->ncol == pX_p->ncol == 9
   
    //A[1][:]* X[:][1]
    pX_m->val[nCol*sigmaIdx + idxSt1] = (A[1][1] * pX_p->val[nCol*idxSt1 + sigmaIdx]) + (A[1][3] * pX_p->val[nCol*idxSt3 + sigmaIdx]);


}

void Fx2(tMatrix * pu_p, tMatrix * pX_p, tMatrix * pX_m,int sigmaIdx)
{
    int nCol = pX_m->ncol; //pX_m->ncol == pX_p->ncol == 9
    
    //fx2() = A[2][:]* X[:][2]
    pX_m->val[nCol*sigmaIdx + idxSt2] = (A[2][2] * pX_p->val[nCol*idxSt2 + sigmaIdx]);
}

void Fx3(tMatrix * pu_p, tMatrix * pX_p, tMatrix * pX_m,int sigmaIdx)
{
    int nCol = pX_m->ncol; //pX_m->ncol == pX_p->ncol == 9
    
    //A[3][:]* X[:][3]
    pX_m->val[nCol*sigmaIdx + idxSt3] = (A[3][3] * pX_p->val[nCol*idxSt3 + sigmaIdx]);
}

void Hy1(tMatrix * pu, tMatrix * pX_m, tMatrix * pY_m,int sigmaIdx)
{
    
}

void Hy2(tMatrix * pu, tMatrix * pX_m, tMatrix * pY_m,int sigmaIdx)
{
    
}
