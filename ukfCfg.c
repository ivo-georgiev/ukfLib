//This file depend of phisical system model

#include "ukfCfg.h"
#include "ukfLib.h"
#include "math.h"

void Fx1(tMatrix * pu_p, tMatrix * px_p, tMatrix * pX_m,int sigmaIdx);
void Fx2(tMatrix * pu_p, tMatrix * px_p, tMatrix * pX_m,int sigmaIdx);
void Fx3(tMatrix * pu_p, tMatrix * px_p, tMatrix * pX_m,int sigmaIdx);
void Fx4(tMatrix * pu_p, tMatrix * px_p, tMatrix * pX_m,int sigmaIdx);

void Hy1(tMatrix * pu, tMatrix * pX_m, tMatrix * pY_m,int sigmaIdx);
void Hy2(tMatrix * pu, tMatrix * pX_m, tMatrix * pY_m,int sigmaIdx);
void Hy3(tMatrix * pu, tMatrix * pX_m, tMatrix * pY_m,int sigmaIdx);
void Hy4(tMatrix * pu, tMatrix * pX_m, tMatrix * pY_m,int sigmaIdx);

static const tPredictFcn PredictFcn[stateVectorLen] = {&Fx1,&Fx2,Fx3,Fx4};
static const tObservFcn  ObservFcn[stateVectorLen] = {&Hy1,&Hy2,Hy3,Hy4};

static const double A[stateVectorLen][stateVectorLen] =
{
    {1,   0, 0.1,   0},
    {0,   1,   0, 0.1},
    {0,   0,   1,   0},
    {0,   0,   0,   1}
};


void Fx1(tMatrix * pu_p, tMatrix * px_p, tMatrix * pX_m,int sigmaIdx)
{
    double * pu = (double *)&pu_p->val;
    int nColXm = pX_m->ncol;
    int nColUp = pu_p->ncol;

    pX_m->val[nColXm*sigmaIdx + idx0] = (A[idx0][idx0] * pu[nColUp*idx0]) + (A[idx0][2] * pu[nColUp * idx2]);

}

void Fx2(tMatrix * pu_p, tMatrix * px_p, tMatrix * pX_m,int sigmaIdx)
{
    double * pu = (double *)&pu_p->val;
    int nColXm = pX_m->ncol;
    int nColUp = pu_p->ncol;
    
    pX_m->val[nColXm*sigmaIdx + idx1] = A[idx1][idx1] * pu[nColUp*idx1] + A[idx1][idx3] * pu[nColUp * idx3];
}

void Fx3(tMatrix * pu_p, tMatrix * px_p, tMatrix * pX_m,int sigmaIdx)
{
    double * pu = (double *)&pu_p->val;
    int nColXm = pX_m->ncol;
    int nColUp = pu_p->ncol;
    
    pX_m->val[nColXm*sigmaIdx + idx2] = A[idx2][idx2] * pu[nColUp * idx2];
}

void Fx4(tMatrix * pu_p, tMatrix * px_p, tMatrix * pX_m,int sigmaIdx)
{
    double * pu = (double *)&pu_p->val;
    int nColXm = pX_m->ncol;
    int nColUp = pu_p->ncol;
    
    pX_m->val[nColXm*sigmaIdx + idx3] = A[idx3][idx3] * pu[nColUp * idx3];
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