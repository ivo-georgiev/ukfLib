#include "System_Types.h"
#include "stdio.h"
#include "math.h"
#include "mtxLib.h"

#ifndef UKFLIB_FILE
#define UKFLIB_FILE

extern const int xMinIdx;
extern const int xMaxIdx;
extern const int xEpsIdx;

extern const int alphaIdx;
extern const int bethaIdx;
extern const int kappaIdx;

typedef void (*PredictFcn_t)(Matrix_t *pu_p, Matrix_t *px_p, Matrix_t *pX_m, int sigmaIdx, double dT);
typedef void (*ObservFcn_t)(Matrix_t *pu, Matrix_t *pX_m, Matrix_t *pY_m, int sigmaIdx);

typedef struct ukfMatrix
{
	Matrix_t Sc_vector;
	Matrix_t Wm_weight_vector;
	Matrix_t Wc_weight_vector;
	Matrix_t x_system_states;
	Matrix_t x_system_states_ic;
	Matrix_t x_system_states_limits;			/**< NOT MANDATORY assign NULL if not required */
	MatrixBool_t x_system_states_limits_enable; /**< NOT MANDATORY assign NULL if not required */
	Matrix_t x_system_states_correction;
	Matrix_t u_system_input;	  /**< NOT MANDATORY assign NULL if not required */
	Matrix_t u_prev_system_input; /**< NOT MANDATORY assign NULL if not required */
	Matrix_t X_sigma_points;
	Matrix_t Y_sigma_points;
	Matrix_t y_predicted_mean;
	Matrix_t y_meas;
	Matrix_t Pyy_out_covariance;
	Matrix_t Pyy_out_covariance_copy;
	Matrix_t Ryy0_init_out_covariance;
	Matrix_t Pxy_cross_covariance;
	Matrix_t Pxx_error_covariance;
	Matrix_t Pxx0_init_error_covariance;
	Matrix_t Qxx_process_noise_cov;
	Matrix_t K_kalman_gain;
#if 0
    Matrix_t K_kalman_gain_transp;
#endif
	Matrix_t I_identity_matrix;
	Matrix_t Pxx_covariance_correction;
	PredictFcn_t *fcnPredict;
	ObservFcn_t *fcnObserve;
	double dT;

} UkfMatrix_t;

typedef struct uKFpar
{
	int xLen; /**< length of state vector */
	int yLen; /**< length of measurement vector */
	int sLen; /**< length of sigma point */
	double alpha; /**< Range:[10e-4 : 1].Smaller alpha leads to a tighter (closer) selection of sigma-points, */
	double betha; /**< Contain information about the prior distribution (for Gaussian, beta = 2 is optimal). */
	double kappa; /**< tertiary scaling parameter, usual value 0. */
	double lambda;
	double dT;
	Matrix_t Wm;
	Matrix_t Wc;
	Matrix_t Qxx;
	Matrix_t Ryy0;
	Matrix_t Pxx0;
	Matrix_t x0;
	Matrix_t xLim;
	MatrixBool_t xLimEnbl;
} UKFpar_t;

typedef struct uKFin
{
	Matrix_t u; /**< \f$ u_k\f$   Current inputs */
	Matrix_t y; /**< \f$ y_k\f$   Current measurement */
} UKFin_t;

typedef struct uKFprev
{
	Matrix_t u_p;   /**< \f$ u_{k-1}\f$   Previous inputs */
	Matrix_t x_p;   /**< \f$ x_{k-1}\f$   Previous states */
	Matrix_t X_p;   /**< \f$ X_{k-1}\f$   Calculate the sigma-points */
	Matrix_t Pxx_p; /**< \f$ P_{k-1}\f$   Previous error covariance  */
} UKFprev_t;

typedef struct uKFpredict /**< p(previous)==k-1, m(minus)=(k|k-1) */
{
	Matrix_t X_m; /**< \f$X_{k|k-1}\f$ Propagate each sigma-point through prediction \f$f(\chi)\f$ */
	Matrix_t x_m; /**< \f$x_{k|k-1}\f$ Calculate mean of predicted state */
	Matrix_t P_m; /**< \f$P_{k|k-1}\f$ Calculate covariance of predicted state   */
	Matrix_t Y_m; /**< \f$Y_{k|k-1}\f$ Propagate each sigma-point through observation */
	Matrix_t y_m; /**< \f$y_{k|k-1}\f$ Calculate mean of predicted output */
	PredictFcn_t *pFcnPredict;
	ObservFcn_t *pFcnObserv;
} UKFpredict_t;

typedef struct uKFupdate
{
	Matrix_t Pyy; /**< Calculate covariance of predicted output */
	Matrix_t Pyy_cpy;
	Matrix_t Pxy; /**< Calculate cross-covariance of state and output */
	Matrix_t K;   /**< \f$K_{k}\f$ Calculate gain */
	Matrix_t x;   /**< \f$x_{k}\f$ Update state estimate */
	Matrix_t x_corr;
	Matrix_t Pxx; /**< \f$P_{k}\f$ Update error covariance */
	Matrix_t Pxx_corr;
	Matrix_t Iyy; /**< tmp buffer initialized as identity matrix stor result from inversion and other operation */
} UKFupdate_t;

typedef struct uKF
{
	UKFpar_t par;
	UKFprev_t prev;
	UKFin_t input;
	UKFpredict_t predict;
	UKFupdate_t update;
} UKF_t;

#endif

extern _Bool ukf_init(UKF_t *const pUkf, UkfMatrix_t *pUkfMatrix);
extern void ukf_step(UKF_t *const pUkf);
