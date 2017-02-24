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
void Hy3(tMatrix * pu, tMatrix * pX_m, tMatrix * pY_m,int sigmaIdx);
void Hy4(tMatrix * pu, tMatrix * pX_m, tMatrix * pY_m,int sigmaIdx);

static const tPredictFcn PredictFcn[stateVectorLen] = {&Fx0,&Fx1,&Fx2,&Fx3};
static const tObservFcn  ObservFcn[stateVectorLen] = {&Hy1,&Hy2,&Hy3,&Hy4};

//-----------------------
//UKF Processing matrix
//-----------------------
double sigma_weight_Wm_1x9[1][9] = {0,0,0,0,0,0,0,0,0};
double sigma_weight_Wc_1x9[1][9] = {0,0,0,0,0,0,0,0,0};

//System input u(k)
double input_u_curr_4x1[4][1] = 
{
    {0},
    {0},
    {0},
    {0}
};
//System input u(k)
double input_u_prev_4x1[4][1] = 
{
    {0},
    {0},
    {0},
    {0}
};
//System states x(k-1)
double state_x_prev_4x1[4][1] = 
{
    {0},
    {0},
    {0},
    {0}
};
//Sigma points X(k-1)
double sigma_X_prev_4x9[4][9]=
{
    {0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0,0},
};
//Sigma points X(k|k-1)
double sigma_X_minus_4x9[4][9]=
{
    {0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0,0},
};

//State covariance X(k|k-1)
double state_cov_P_minus_4x4[4][4]=
{
    {0,0,0,0},
    {0,0,0,0},
    {0,0,0,0},
    {0,0,0,0},
};

//temporal result ,atrix 4*4
double temp_4x4[4][4]=
{
    {0,0,0,0},
    {0,0,0,0},
    {0,0,0,0},
    {0,0,0,0},
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

void Hy3(tMatrix * pu, tMatrix * pX_m, tMatrix * pY_m,int sigmaIdx)
{
    
}

void Hy4(tMatrix * pu, tMatrix * pX_m, tMatrix * pY_m,int sigmaIdx)
{
    
}