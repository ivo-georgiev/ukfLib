#include "ukfLib.h"

extern void Fx0(tMatrix * pu_p, tMatrix * pX_p, tMatrix * pX_m,uint8 sigmaIdx);
extern void Fx1(tMatrix * pu_p, tMatrix * pX_p, tMatrix * pX_m,uint8 sigmaIdx);
extern void Fx2(tMatrix * pu_p, tMatrix * pX_p, tMatrix * pX_m,uint8 sigmaIdx);
extern void Fx3(tMatrix * pu_p, tMatrix * pX_p, tMatrix * pX_m,uint8 sigmaIdx);

extern void Hy1(tMatrix * pu, tMatrix * pX_m, tMatrix * pY_m,uint8 sigmaIdx);
extern void Hy2(tMatrix * pu, tMatrix * pX_m, tMatrix * pY_m,uint8 sigmaIdx);

#define stateVectorLen (uint8)4
#define measVectorLen (uint8)2

#define dT0 (float64)0.1 // [sec] 

extern float64 Wm_sigma_weight_1x9[1][9];
extern float64 Wc_sigma_weight_1x9[1][9];
extern float64 u_curr_system_input_4x1[4][1]; 

//System input previous u(k)
extern float64 u_prev_system_input_4x1[4][1];

//System output measurement y(k)
extern float64 y_curr_system_meas_2x1[2][1];

//System output predicted y_m(k|k-1)
extern float64 y_mean_system_predict_2x1[2][1];

//System states: x(k), x(k-1), x(k|k-1) common array for all 
extern float64 x_system_states_4x1[4][1];

extern float64 x_system_states_ic_4x1[4][1];

extern float64 x_system_states_correction_4x1[4][1];

extern float64 X_sigma_points_4x9[4][9];

extern float64 Y_sigma_points_2x9[2][9];

extern float64 Px_state_cov_4x4[4][4];

extern float64 Pxx_covariance_correction_4x4[4][4];

extern float64 Pyy_out_cov_2x2[2][2];

extern float64 Pyy_out_cov_copy_2x2[2][2];

extern float64 Pxy_state_out_cov_4x2[4][2];

extern float64 K_kalman_gain_4x2[4][2];

extern float64 K_kalman_transp_gain_2x4[2][4];

extern float64 Qxx_process_noise_cov_4x4[4][4];

extern float64 temporal_2x2[2][2];

extern float64 P0_state_cov_4x4[4][4];

extern float64 Ryy_out_cov_noise_2x2[2][2];

extern tPredictFcn PredictFcn[stateVectorLen];

extern tObservFcn  ObservFcn[measVectorLen];