#include "mtxLib.h"

#ifndef UKFLIB_FILE
#define UKFLIB_FILE

typedef void (* tPredictFcn) (tMatrix * pu_p, tMatrix * px_p, tMatrix * pX_m,int sigmaIdx);
typedef void (* tObservFcn) (tMatrix * pu, tMatrix * pX_m, tMatrix * pY_m,int sigmaIdx);

typedef struct tUkfMatrix
{
    tMatrix Wm_weight_vector;
    tMatrix Wc_weight_vector;
    tMatrix x_system_states;
    tMatrix x_system_states_ic;
    tMatrix x_system_states_correction;
    tMatrix u_system_input;
    tMatrix u_prev_system_input;
    tMatrix X_sigma_points;
    tMatrix Y_sigma_points;
    tMatrix y_predicted_mean;
    tMatrix y_meas;
    tMatrix Pyy_out_covariance;
    tMatrix Pyy_out_covariance_copy;
    tMatrix Ryy0_init_out_covariance;
    tMatrix Pxy_cross_covariance;
    tMatrix Pxx_error_covariance;
    tMatrix Pxx0_init_error_covariance;
    tMatrix Qxx_process_noise_cov;
    tMatrix K_kalman_gain;
    tMatrix K_kalman_gain_transp;
    tMatrix I_identity_matrix;
    tMatrix Pxx_covariance_correction; 
    tPredictFcn * fcnPredict;
    tObservFcn * fcnObserve;

}tUkfMatrix;



#define alphaIdx   (int)0
#define bethaIdx   (int)1
#define kappaIdx   (int)2
#define scalingLen  (int)3


typedef struct tUKFpar
{
    int xLen;//length of state vector
    int yLen;//length of measurement vector
    int sLen;//length of sigma point
    double alpha;//Range:[10e-4 : 1].Smaller alpha leads to a tighter (closer) selection of sigma-points,
    double betha;//Contain information about the prior distribution (for Gaussian, beta = 2 is optimal).
    double kappa; //tertiary scaling parameter, usual value 0.
    double lambda;
    tMatrix Wm;
    tMatrix Wc;
    tMatrix Qxx;
    tMatrix Ryy0;
    tMatrix Pxx0;
    tMatrix x0;

}tUKFpar;

typedef struct tUKFin
{
    tMatrix u;    // u(k)   Current inputs
    tMatrix y;    // y(k)   Current measurement
}tUKFin;

typedef struct tUKFprev
{
    tMatrix u_p;    // u(k-1)   Previous inputs
    tMatrix x_p;    // x(k-1)   Previous states
    tMatrix X_p;    // X(k-1)   Calculate the sigma-points
    tMatrix Pxx_p;    // P(k-1)    Previous error covariance 
}tUKFprev;

typedef struct tUKFpredict //p(previous)==k-1, m(minus)=(k|k-1)
{
    tMatrix X_m;    //X(k|k-1) Propagate each sigma-point through prediction f(Chi)
    tMatrix x_m;    //x(k|k-1) Calculate mean of predicted state
    tMatrix P_m;    //P(k|k-1) Calculate covariance of predicted state  
    tMatrix Y_m;    //Y(k|k-1) Propagate each sigma-point through observation
    tMatrix y_m;    //y(k|k-1) Calculate mean of predicted output
    tPredictFcn * pFcnPredict;
    tObservFcn * pFcnObserv;
}tUKFpredict;



typedef struct tUKFupdate
{
    tMatrix Pyy;    //Calculate covariance of predicted output
    tMatrix Pyy_cpy; 
    tMatrix Pxy;   //Calculate cross-covariance of state and output
    tMatrix K;     //K(k) Calculate gain
    tMatrix Kt;     //Kt(k) Kalman gain transponce
    tMatrix x;     //x(k) Update state estimate
    tMatrix x_corr;
    tMatrix Pxx;     //P(k) Update error covariance
    tMatrix Pxx_corr;
    tMatrix Iyy;     //tmp buffer initialized as identity matrix stor result from inversion and other operation  
}tUKFupdate;


typedef struct tUKF
{
    tUKFpar     par;
    tUKFprev    prev;
    tUKFin      input;
    tUKFpredict predict;
    tUKFupdate  update;

}tUKF;

#endif

extern int ukf_init(tUKF * const pUkf,double scaling[scalingLen],int xLen,int yLen, tUkfMatrix * pUkfMatrix);
extern void ukf_step(tUKF * const pUkf);