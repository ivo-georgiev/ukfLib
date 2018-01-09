 /******************************************************************************************************************************************************************************************************\
 *** 
 *** Description       : IMPLEMENTATION OF THE ADDITIVE NOISE UKF: 
 *** Codefile          : ukfCfgTemplate.c
 *** Documentation     : 
 ***
 *** MIT License
 ***
 *** Copyright (c) 2018 ivo-georgiev
 ***  
 *** Permission is hereby granted, free of charge, to any person obtaining a copy
 *** of this software and associated documentation files (the "Software"), to deal
 *** in the Software without restriction, including without limitation the rights
 *** to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 *** copies of the Software, and to permit persons to whom the Software is
 *** furnished to do so, subject to the following conditions:
 ***    
 *** The above copyright notice and this permission notice shall be included in all
 *** copies or substantial portions of the Software.
 ***      
 *** THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 *** IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 *** FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 *** AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 *** LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 *** OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 *** SOFTWARE.
\******************************************************************************************************************************************************************************************************/
#include "ukfCfg<cfgId>.h"

//<STATE TRANSITION PROTOTYPE:BEGIN>
//<STATE TRANSITION PROTOTYPE:END>

//<MEASUREMENT PROTOTYPE:BEGIN>
//<MEASUREMENT PROTOTYPE:END>

//<STATE TRANSITION PTR ARRAY:BEGIN>
//<STATE TRANSITION PTR ARRAY:END>

//<MEASUREMENT PTR ARRAY:BEGIN>
//<MEASUREMENT PTR ARRAY:END>

//-----------------------
//UKF Processing matrix
//-----------------------
static float64 Sc_vector[1][3] = {{alpha,betha,kappa}};
static float64 Wm_weight_vector[1][sL] = {{<>}};
static float64 Wc_weight_vector[1][sL] = {{<>}};
//static float64 u_system_input[2][1] = {{0},{0}}; 
//static float64 u_prev_system_input[2][1] = {{0},{0}};
static float64 y_meas[yL][1] = <matrix> ;
static float64 y_predicted_mean[yL][1] = <matrix> ;
static float64 x_system_states[xL][1] = <matrix>;
static float64 x_system_states_ic[xL][1] = <matrix>;
static float64 x_system_states_limits[xL][3] = <matrix>;
static boolean x_system_states_limits_enable[xL][1] = <matrix>;
static float64 x_system_states_correction[xL][1] = <matrix>;
static float64 X_sigma_points[xL][sL]= <matrix>;
static float64 Y_sigma_points[yL][sL]= <matrix>;
static float64 Pxx_error_covariance[xL][xL]= <matrix>;
static float64 Pxx0_init_error_covariance[xL][xL]= <matrix>;
static float64 Qxx_process_noise_cov[xL][xL]= <matrix>;
static float64 Ryy0_init_out_covariance[yL][yL]= <matrix>;
static float64 Pyy_out_covariance[yL][yL]= <matrix>;
static float64 Pyy_out_covariance_copy[yL][yL]= <matrix>;
static float64 Pxy_cross_covariance[xL][yL]= <matrix>;
static float64 K_kalman_gain[xL][yL]= <matrix>;
static float64 Pxx_covariance_correction[xL][yL]= <matrix>;
static float64 I_identity_matrix[yL][yL]= <matrix>;

tUkfMatrix UkfMatrixCfg<cfgId> = 
{
    {COLXROW(Sc_vector),NROWS(Sc_vector),NCOL(Sc_vector),&Sc_vector[0][0]},
    {COLXROW(Wm_weight_vector),NROWS(Wm_weight_vector),NCOL(Wm_weight_vector),&Wm_weight_vector[0][0]},
    {COLXROW(Wc_weight_vector),NROWS(Wc_weight_vector),NCOL(Wc_weight_vector),&Wc_weight_vector[0][0]},
    {COLXROW(x_system_states),NROWS(x_system_states),NCOL(x_system_states),&x_system_states[0][0]},
    {COLXROW(x_system_states_ic),NROWS(x_system_states_ic),NCOL(x_system_states_ic),&x_system_states_ic[0][0]},
    {COLXROW(x_system_states_limits),NROWS(x_system_states_limits),NCOL(x_system_states_limits),&x_system_states_limits[0][0]},
    {COLXROW(x_system_states_limits_enable),NROWS(x_system_states_limits_enable),NCOL(x_system_states_limits_enable),&x_system_states_limits_enable[0][0]},
    {COLXROW(x_system_states_correction),NROWS(x_system_states_correction),NCOL(x_system_states_correction),&x_system_states_correction[0][0]},
    {0,0,0,NULL},//{COLXROW(u_system_input),NROWS(u_system_input),NCOL(u_system_input),&u_system_input[0][0]},
    {0,0,0,NULL},//{COLXROW(u_prev_system_input),NROWS(u_prev_system_input),NCOL(u_prev_system_input),&u_prev_system_input[0][0]},
    {COLXROW(X_sigma_points),NROWS(X_sigma_points),NCOL(X_sigma_points),&X_sigma_points[0][0]},
    {COLXROW(Y_sigma_points),NROWS(Y_sigma_points),NCOL(Y_sigma_points),&Y_sigma_points[0][0]},
    {COLXROW(y_predicted_mean),NROWS(y_predicted_mean),NCOL(y_predicted_mean),&y_predicted_mean[0][0]},   
    {COLXROW(y_meas),NROWS(y_meas),NCOL(y_meas),&y_meas[0][0]},
    {COLXROW(Pyy_out_covariance),NROWS(Pyy_out_covariance),NCOL(Pyy_out_covariance),&Pyy_out_covariance[0][0]},
    {COLXROW(Pyy_out_covariance_copy),NROWS(Pyy_out_covariance_copy),NCOL(Pyy_out_covariance_copy),&Pyy_out_covariance_copy[0][0]},
    {COLXROW(Ryy0_init_out_covariance),NROWS(Ryy0_init_out_covariance),NCOL(Ryy0_init_out_covariance),&Ryy0_init_out_covariance[0][0]},
    {COLXROW(Pxy_cross_covariance),NROWS(Pxy_cross_covariance),NCOL(Pxy_cross_covariance),&Pxy_cross_covariance[0][0]},
    {COLXROW(Pxx_error_covariance),NROWS(Pxx_error_covariance),NCOL(Pxx_error_covariance),&Pxx_error_covariance[0][0]},
    {COLXROW(Pxx0_init_error_covariance),NROWS(Pxx0_init_error_covariance),NCOL(Pxx0_init_error_covariance),&Pxx0_init_error_covariance[0][0]},
    {COLXROW(Qxx_process_noise_cov),NROWS(Qxx_process_noise_cov),NCOL(Qxx_process_noise_cov),&Qxx_process_noise_cov[0][0]},
    {COLXROW(K_kalman_gain),NROWS(K_kalman_gain),NCOL(K_kalman_gain),&K_kalman_gain[0][0]},
    {COLXROW(I_identity_matrix),NROWS(I_identity_matrix),NCOL(I_identity_matrix),&I_identity_matrix[0][0]},  
    {COLXROW(Pxx_covariance_correction),NROWS(Pxx_covariance_correction),NCOL(Pxx_covariance_correction),&Pxx_covariance_correction[0][0]},   
    &PredictFcn[0],
    &ObservFcn[0],
    dT
};

//<STATE TRANSITION:BEGIN>
//<STATE TRANSITION:END>

//<MEASUREMENT FUNCTION:BEGIN>
//<MEASUREMENT FUNCTION:END>

