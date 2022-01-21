#include "System_Types.h"
#include "stdio.h"
#include "math.h"
#include "mtxLib.h"

#ifndef UKFLIB_FILE
#define UKFLIB_FILE

#define xMinIdx (uint8_t)0
#define xMaxIdx (uint8_t)1
#define xEpsIdx (uint8_t)2

#define alphaIdx   (uint8_t)0
#define bethaIdx   (uint8_t)1
#define kappaIdx   (uint8_t)2

typedef void (* tPredictFcn) (tMatrix * pu_p, tMatrix * px_p, tMatrix * pX_m,uint8_t sigmaIdx, float64 dT);
typedef void (* tObservFcn) (tMatrix * pu, tMatrix * pX_m, tMatrix * pY_m,uint8_t sigmaIdx);

typedef struct ukfMatrix
{
    tMatrix Sc_vector;
    tMatrix Wm_weight_vector;
    tMatrix Wc_weight_vector;
    tMatrix x_system_states;
    tMatrix x_system_states_ic;
    tMatrix x_system_states_limits;              //NOT MANDATORY assign NULL if not required
    tMatrixBool x_system_states_limits_enable;   //NOT MANDATORY assign NULL if not required
    tMatrix x_system_states_correction;
    tMatrix u_system_input;                      //NOT MANDATORY assign NULL if not required
    tMatrix u_prev_system_input;                 //NOT MANDATORY assign NULL if not required
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
    //tMatrix K_kalman_gain_transp;
    tMatrix I_identity_matrix;
    tMatrix Pxx_covariance_correction; 
    tPredictFcn * fcnPredict;
    tObservFcn * fcnObserve;
    float64 dT;

}tUkfMatrix;

typedef struct uKFpar
{
    uint8_t xLen;//length of state vector
    uint8_t yLen;//length of measurement vector
    uint8_t sLen;//length of sigma point
    float64 alpha;//Range:[10e-4 : 1].Smaller alpha leads to a tighter (closer) selection of sigma-points,
    float64 betha;//Contain information about the prior distribution (for Gaussian, beta = 2 is optimal).
    float64 kappa; //tertiary scaling parameter, usual value 0.
    float64 lambda;
    float64 dT;
    tMatrix Wm;
    tMatrix Wc;
    tMatrix Qxx;
    tMatrix Ryy0;
    tMatrix Pxx0;
    tMatrix x0;
    tMatrix xLim;
    tMatrixBool xLimEnbl;
}tUKFpar;

typedef struct uKFin
{
    tMatrix u;    // u(k)   Current inputs
    tMatrix y;    // y(k)   Current measurement
}tUKFin;

typedef struct uKFprev
{
    tMatrix u_p;    // u(k-1)   Previous inputs
    tMatrix x_p;    // x(k-1)   Previous states
    tMatrix X_p;    // X(k-1)   Calculate the sigma-points
    tMatrix Pxx_p;  // P(k-1)    Previous error covariance 
}tUKFprev;

typedef struct uKFpredict //p(previous)==k-1, m(minus)=(k|k-1)
{
    tMatrix X_m;    //X(k|k-1) Propagate each sigma-point through prediction f(Chi)
    tMatrix x_m;    //x(k|k-1) Calculate mean of predicted state
    tMatrix P_m;    //P(k|k-1) Calculate covariance of predicted state  
    tMatrix Y_m;    //Y(k|k-1) Propagate each sigma-point through observation
    tMatrix y_m;    //y(k|k-1) Calculate mean of predicted output
    tPredictFcn * pFcnPredict;
    tObservFcn * pFcnObserv;
}tUKFpredict;

typedef struct uKFupdate
{
    tMatrix Pyy;    //Calculate covariance of predicted output
    tMatrix Pyy_cpy; 
    tMatrix Pxy;     //Calculate cross-covariance of state and output
    tMatrix K;       //K(k) Calculate gain
    tMatrix x;       //x(k) Update state estimate
    tMatrix x_corr;
    tMatrix Pxx;     //P(k) Update error covariance
    tMatrix Pxx_corr;
    tMatrix Iyy;     //tmp buffer initialized as identity matrix stor result from inversion and other operation  
}tUKFupdate;

typedef struct uKF
{
    tUKFpar     par;
    tUKFprev    prev;
    tUKFin      input;
    tUKFpredict predict;
    tUKFupdate  update;
}tUKF;

#endif

extern _Bool ukf_init(tUKF * const pUkf, tUkfMatrix * pUkfMatrix);
extern void ukf_step(tUKF * const pUkf);
