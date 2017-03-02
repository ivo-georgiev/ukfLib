
#define idxSt0 (int)0
#define idxSt1 (int)1
#define idxSt2 (int)2
#define idxSt3 (int)3

#define stateVectorLen (int)4

#define dT0 (double)0.1 // [sec] 

extern double Wm_sigma_weight_1x9[9];
extern double Wc_sigma_weight_1x9[9];
extern double u_curr_system_input_4x1[4]; 

//System input previous u(k)
extern double u_prev_system_input_4x1[4];

//System output measurement y(k)
extern double y_curr_system_meas_2x1[4];

//System output predicted y_m(k|k-1)
extern double y_mean_system_predict_2x1[4];

//System states: x(k), x(k-1), x(k|k-1) common array for all 
extern double x_system_states_4x1[4];

extern double X_sigma_points_4x9[4][9];

extern double Y_sigma_points_2x9[2][9];

extern double Px_state_cov_4x4[4][4];

extern double Pyy_out_cov_2x2[2][2];


extern double Pxy_state_out_cov_4x2[4][2];

extern double K_kalman_gain_4x2[4][2];

extern double K_kalman_transp_gain_4x2[2][4];

extern double temporal_4x4[4][4];