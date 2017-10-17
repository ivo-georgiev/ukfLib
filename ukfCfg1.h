#ifndef UKFCFG1_H
#define UKFCFG1_H

#include "ukfLib.h"

//#define dT0 (float64)0.1 // [sec] 

extern float64 _Sc_vector_1x3[1][3];
extern float64 _Wm_sigma_weight_1x9[1][9];
extern float64 _Wc_sigma_weight_1x9[1][9];
extern float64 _u_curr_system_input_4x1[4][1]; 
extern float64 _u_prev_system_input_4x1[4][1];
extern float64 _y_curr_system_meas_2x1[2][1];
extern float64 _y_mean_system_predict_2x1[2][1]; 
extern float64 _x_system_states_4x1[4][1];
extern float64 _x_system_states_ic_4x1[4][1];
extern float64 _x_system_states_correction_4x1[4][1];
extern float64 _X_sigma_points_4x9[4][9];
extern float64 _Y_sigma_points_2x9[2][9];
extern float64 _Px_state_cov_4x4[4][4];
extern float64 _Pxx_covariance_correction_4x4[4][4];
extern float64 _Pyy_out_cov_2x2[2][2];
extern float64 _Pyy_out_cov_copy_2x2[2][2];
extern float64 _Pxy_state_out_cov_4x2[4][2];
extern float64 _K_kalman_gain_4x2[4][2];
extern float64 _K_kalman_transp_gain_2x4[2][4];
extern float64 _Qxx_process_noise_cov_4x4[4][4];
extern float64 _temporal_2x2[2][2];
extern float64 _P0_state_cov_4x4[4][4];
extern float64 _Ryy_out_cov_noise_2x2[2][2];

extern tPredictFcn _PredictFcn[4];

extern tObservFcn  _ObservFcn[2];
#endif /* UKFCFG_H */

