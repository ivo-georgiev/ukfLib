#include "System_Types.h"
#include "stdio.h"
#include "math.h"
#include "mtxLib.h"

#ifndef UKFLIB_FILE
#define UKFLIB_FILE

extern const uint8_t xMinIdx;
extern const uint8_t xMaxIdx;
extern const uint8_t xEpsIdx;

extern const uint8_t alphaIdx;
extern const uint8_t bethaIdx;
extern const uint8_t kappaIdx;

typedef void (* PredictFcn_t) (Matrix64_t * pu_p, Matrix64_t * px_p, Matrix64_t * pX_m,uint8_t sigmaIdx, float64 dT);
typedef void (* ObservFcn_t) (Matrix64_t * pu, Matrix64_t * pX_m, Matrix64_t * pY_m,uint8_t sigmaIdx);

typedef struct ukfMatrix
{
    Matrix64_t Sc_vector;
    Matrix64_t Wm_weight_vector;
    Matrix64_t Wc_weight_vector;
    Matrix64_t x_system_states;
    Matrix64_t x_system_states_ic;
    Matrix64_t x_system_states_limits;              /**< NOT MANDATORY assign NULL if not required */
    MatrixBool64_t x_system_states_limits_enable;   /**< NOT MANDATORY assign NULL if not required */
    Matrix64_t x_system_states_correction;
    Matrix64_t u_system_input;                      /**< NOT MANDATORY assign NULL if not required */
    Matrix64_t u_prev_system_input;                 /**< NOT MANDATORY assign NULL if not required */
    Matrix64_t X_sigma_points;
    Matrix64_t Y_sigma_points;
    Matrix64_t y_predicted_mean;
    Matrix64_t y_meas;
    Matrix64_t Pyy_out_covariance;
    Matrix64_t Pyy_out_covariance_copy;
    Matrix64_t Ryy0_init_out_covariance;
    Matrix64_t Pxy_cross_covariance;
    Matrix64_t Pxx_error_covariance;
    Matrix64_t Pxx0_init_error_covariance;
    Matrix64_t Qxx_process_noise_cov;
    Matrix64_t K_kalman_gain;
#if 0
    Matrix64_t K_kalman_gain_transp;
#endif
    Matrix64_t I_identity_matrix;
    Matrix64_t Pxx_covariance_correction;
    PredictFcn_t * fcnPredict;
    ObservFcn_t * fcnObserve;
    float64 dT;

}UkfMatrix64_t;

typedef struct uKFpar
{
    uint8_t xLen;/**< length of state vector */
    uint8_t yLen;/**< length of measurement vector */
    uint8_t sLen;/**< length of sigma point */
    float64 alpha;/**< Range:[10e-4 : 1].Smaller alpha leads to a tighter (closer) selection of sigma-points, */
    float64 betha;/**< Contain information about the prior distribution (for Gaussian, beta = 2 is optimal). */
    float64 kappa; /**< tertiary scaling parameter, usual value 0. */
    float64 lambda;
    float64 dT;
    Matrix64_t Wm;
    Matrix64_t Wc;
    Matrix64_t Qxx;
    Matrix64_t Ryy0;
    Matrix64_t Pxx0;
    Matrix64_t x0;
    Matrix64_t xLim;
    MatrixBool64_t xLimEnbl;
}UKFpar64_t;

typedef struct uKFin
{
    Matrix64_t u;    /**< \f$ u_k\f$   Current inputs */
    Matrix64_t y;    /**< \f$ y_k\f$   Current measurement */
}UKFin64_t;

typedef struct uKFprev
{
    Matrix64_t u_p;    /**< \f$ u_{k-1}\f$   Previous inputs */
    Matrix64_t x_p;    /**< \f$ x_{k-1}\f$   Previous states */
    Matrix64_t X_p;    /**< \f$ X_{k-1}\f$   Calculate the sigma-points */
    Matrix64_t Pxx_p;  /**< \f$ P_{k-1}\f$   Previous error covariance  */
}UKFprev64_t;

typedef struct uKFpredict /**< p(previous)==k-1, m(minus)=(k|k-1) */
{
    Matrix64_t X_m;    /**< \f$X_{k|k-1}\f$ Propagate each sigma-point through prediction \f$f(\chi)\f$ */
    Matrix64_t x_m;    /**< \f$x_{k|k-1}\f$ Calculate mean of predicted state */
    Matrix64_t P_m;    /**< \f$P_{k|k-1}\f$ Calculate covariance of predicted state   */
    Matrix64_t Y_m;    /**< \f$Y_{k|k-1}\f$ Propagate each sigma-point through observation */
    Matrix64_t y_m;    /**< \f$y_{k|k-1}\f$ Calculate mean of predicted output */
    PredictFcn_t * pFcnPredict;
    ObservFcn_t * pFcnObserv;
}UKFpredict64_t;

typedef struct uKFupdate
{
    Matrix64_t Pyy;    /**< Calculate covariance of predicted output */
    Matrix64_t Pyy_cpy;
    Matrix64_t Pxy;     /**< Calculate cross-covariance of state and output */
    Matrix64_t K;       /**< \f$K_{k}\f$ Calculate gain */
    Matrix64_t x;       /**< \f$x_{k}\f$ Update state estimate */
    Matrix64_t x_corr;
    Matrix64_t Pxx;     /**< \f$P_{k}\f$ Update error covariance */
    Matrix64_t Pxx_corr;
    Matrix64_t Iyy;     /**< tmp buffer initialized as identity matrix stor result from inversion and other operation   */
}UKFupdate64_t;

typedef struct uKF
{
    UKFpar64_t     par;
    UKFprev64_t    prev;
    UKFin64_t      input;
    UKFpredict64_t predict;
    UKFupdate64_t  update;
}UKF64_t;

#endif

extern _Bool ukf_init(UKF64_t * const pUkf, UkfMatrix64_t * pUkfMatrix);
extern void ukf_step(UKF64_t * const pUkf);
