#include "ukfCfg.h"

/*---------------------------------------------*/
/*           Function Prototype                */
/*---------------------------------------------*/
void show_matrix_obj(tMatrix A);
void show_matrix(float64 * A, int n,int m);
void ukf_test(void);

//


/*---------------------------------------------*/
/*             Sample Matrix                   */
/*---------------------------------------------*/
float64 symMtx[5][5] =
{{ 3.009964742796415,  -1.009481719676593,   1.161774056277037,  -0.291928007896218,  -0.775215823251770},
{ -1.009481719676593,   3.009964742796415,  -1.009481719676593,   1.161774056277037,  -0.291928007896218},
{  1.161774056277037,  -1.009481719676593,   3.009964742796415,  -1.009481719676593,   1.161774056277037},
{ -0.291928007896218,   1.161774056277037,  -1.009481719676593,   3.009964742796415,  -1.009481719676593},
{ -0.775215823251770,  -0.291928007896218,   1.161774056277037,  -1.009481719676593,   3.009964742796415}};

float64 symMtxChol[5][5] =
{{1.734924996302842,     -0.581858997840148,   0.669639355449256,      -0.168265491890613,  -0.446829589119858},
{0                 ,      1.634443284249678,  -0.379239855780691,       0.650904975441148,  -0.337680621986338},
{0                 ,      0                ,   1.554903536627710,      -0.418003689501540,   0.857240820764834},
{0                 ,      0                ,   0                ,       1.543776893059448,  -0.328117294491480},
{0                 ,      0                ,   0                ,       0                ,   1.361527478565284}};

float64 Identity_5x5[4][4] =
{{1.0,  0     ,0,     0,     },
{0,     1.0   ,0,     0,     },
{0,     0     ,1.0,   0,     },
{0,     0     ,0,     1.0   }};

float64 TestMatrix_0_4x4[4][4] =
{{3.0,  5.0, -1.0,  -4},
{ 1.0,  4.0, -0.7,  -3},
{ 0,   -2.0,  0,     1},
{-2.0,  6.0,  0,     0.3}};


float64 TestMatrix_1_3x3[3][3] = 
{{10.5, 2.17, 3.03},
{ 0.44, 0.59, 6.89},
{ 7.56, 8.17, 9.21}};

float64 TestMatrix_2_3x3[3][3] = 
{ {1.11, 29.3, 31.2},
{45.3, 5.17, 6.11},
{7.61, 88.0, 9.34}};


float64 TestMatrix_1_2x3[2][3] = 
{ {1.11, 29.3, 31.2},   //size 3x3
  {45.3, 5.17, 6.11}};

float64 TestMatrixDest_3x2[3][2] =
{
    {0,0},
    {0,0},
    {0,0}  
};


void main(void)
{
    tMatrix myFactMatrix;
    tMatrix myTestMatx={0,0,NULL};
    tMatrix myChol={0,0,NULL};
    tMatrix Im={0,0,NULL};
    tMatrix oTestMatrix_0_4x4={0,0,NULL};
    tMatrix oTestMatrixDest_3x2={0,0,NULL};

    mtx_init_f64(&myTestMatx,&TestMatrix_1_2x3[0][0],NROWS(TestMatrix_1_2x3),NCOL(TestMatrix_1_2x3));
    mtx_init_f64(&myChol,&symMtxChol[0][0],NROWS(symMtxChol),NCOL(symMtxChol));
    mtx_init_f64(&Im,&Identity_5x5[0][0],NROWS(Identity_5x5),NCOL(Identity_5x5));
    mtx_init_f64(&oTestMatrix_0_4x4,&TestMatrix_0_4x4[0][0],NROWS(TestMatrix_0_4x4),NCOL(TestMatrix_0_4x4));
    mtx_init_f64(&oTestMatrixDest_3x2,&TestMatrixDest_3x2[0][0],NROWS(TestMatrixDest_3x2),NCOL(TestMatrixDest_3x2));
    //show_matrix(&TestMatrix_1_2x3[0][0],2,3);
    //show_matrix_obj(myTestMatx);

    (void)mtx_init_f64(&myFactMatrix,&symMtx[0][0],NROWS(symMtx),NCOL(symMtx));
    //show_matrix(&symMtx[0][0],5,5);

    (void)mtx_chol_f64(&myFactMatrix);
    //show_matrix_obj(myFactMatrix);

    show_matrix_obj(myChol);
    mtx_transp_square_f64(&myChol);
    show_matrix_obj(myChol);


    show_matrix_obj(Im);
    show_matrix_obj(oTestMatrix_0_4x4);
    mtx_inv_f64(&oTestMatrix_0_4x4, &Im);
    show_matrix_obj(Im);

    show_matrix_obj(myTestMatx);
    mtx_transp_dest_f64(&myTestMatx,&oTestMatrixDest_3x2);
    show_matrix_obj(oTestMatrixDest_3x2);

    //UKF test start here
    ukf_test();

}

void show_matrix_obj(tMatrix A)
{
    int i,j;
    
    for (i = 0; i < A.nrow; i++) 
    {
        for (j = 0; j < A.ncol; j++)
        {
            printf("%2.5f ", A.val[A.ncol*i + j]);
        }
        printf("\n");
    }
    printf("\n");
}

void show_matrix(float64 * A, int n,int m)
{
    int i,j;
    
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < m; j++)
        {
            printf("%2.5f ", A[m*i + j]);
        }
        printf("\n");
    }
    printf("\n");
}
/******************************************************************************************************************************************************************************************************\
 ***  FUNCTION:
 ***      void ukf_test(void)
 *** 
 ***  DESCRIPTION:
 ***      Initialize and test UKF C implementation against expected result. Filter is tested in the loop from 15 steps. 
 ***      Total root square error is accumulated in the same loop for each state in order to show deviation from reference matlab solution.      
 ***            
 ***  PARAMETERS:
 ***      Type               Name              Range              Description
 ***      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 ***      void
 ***  RETURNS:
 ***      void
 ***
 ***  SETTINGS:
 ***
\******************************************************************************************************************************************************************************************************/
void ukf_test(void)
{
    boolean tfInitFail = 0;
    float64 param[3]={1,2,0};
    tUkfMatrix UkfMat;
    tUKF ukfIo;
    tMatrix myTestMatx;
    uint8 simLoop;

    static const float64 yt[2][15]=
    {
        {0,  16.085992708563385,  12.714829185978214,  14.528500994457660,  19.105561355310275, 23.252820029388918,  29.282949862903255,  36.270058819651275,  44.244884173240955,  47.394243121124411, 55.988905459180458, 61.667450941562109, 68.624980301613647, 76.337963872393104, 82.611325690835159}, //y1 test
        {0,  16.750821420874981,  14.277640835870006,  16.320754051600520,  20.560460303503849, 24.827446289454556,  31.290961393448615,  36.853553457560210,  42.157283183453522,  49.382835230961490, 57.516319669684677, 65.664496283509095, 71.428712755732704, 79.241720894223079, 84.902760328915676 }
    };

    static const float64 xt[4][15]=
    {
        {-0.533020535100940, 4.504298532456061,   9.695289204302437, 14.807408774583198, 19.790259984337496, 24.742668128680663,  29.820581728530414, 34.967601515224246, 40.413443563097616, 46.284515939371857, 52.277685414212016,  58.085505950486329,  63.608578996636339,  69.186393953977429,  75.109287318587334},
        {-0.434369839833425, 4.365037434994960,   9.205660222476737, 14.025603873653782, 18.783261077812174, 23.384002685293172,  28.136582816181591, 32.599696482981962, 36.771914052720703, 41.250093329804031, 45.608894967716552,  49.989876036228665,  54.220796522040331,  58.377142592484368,  62.501946864068032},
        {50.373190675570008, 51.909906718463759, 51.121195702807604, 49.828512097542969, 49.524081443431683, 50.779135998497523,  51.470197866938307, 54.458420478733728, 58.710723762742383, 59.931694748401569, 58.078205362743148,  55.230730461500130,  55.778149573410907,  59.228933646099108,  61.382169390093580},
        {47.994072748283855, 48.406227874817773, 48.199436511770458, 47.576572041583923, 46.007416074809989, 47.525801308884184,  44.631136668003741, 41.722175697387378, 44.781792770833249, 43.588016379125193, 43.809810685121128,  42.309204858116686,  41.563460704440388,  41.248042715836618,  43.019251841426751}

    };

    static const float64 x_exp[15][4]=
    { /*          x1                   x2                   x3                   x4*/
        {4.901482729572258,    4.576939885855807,  49.990342921246459,  49.958134463327802},
        {10.103304943868373,   9.409135720815829,  50.226544716205318,  49.750795004242228},
        {15.132069573131298,  14.138974122835807,  50.429540890147599,  49.191128327737864},
        {20.322823824348411,  19.096919763991380,  50.836860772439010,  49.189580207886742},
        {24.940146120267713,  23.399758647105461,  49.577386595072561,  47.383813382660449},
        {30.021901202161214,  27.882145089050120,  49.977568320123794,  46.551744562547626},
        {34.844137036519108,  32.753891693435087,  49.474006205027358,  47.190693993547214},
        {39.048783329419251,  38.499098203031146,  47.606199725375902,  50.001113730363919},
        {43.883085498256158,  42.383331307689538,  47.209657232695072,  46.747611757031784},
        {49.479941190207498,  47.255980559687778,  49.911944505272395,  47.887236233284476},
        {55.928745858553086,  51.180270357916882,  53.472542964944132,  46.510558543249353},
        {61.636426955126616,  55.275415649334157,  54.052126522632797,  45.262815392265203},
        {67.755622369652016,  59.602096868661732,  55.881393486598796,  45.326766104509289},
        {73.045763444967164,  63.838187852739992,  54.782159791340007,  44.291415099856643},
        {80.489525793047093,  66.908477563332085,  58.973616985147245,  42.638148924845950}
    };

    float64 rootSquareErr_X0 = 0;
    float64 rootSquareErr_X1 = 0;
    float64 rootSquareErr_X2 = 0;
    float64 rootSquareErr_X3 = 0;



    mtx_init_f64(&myTestMatx,&TestMatrix_1_2x3[0][0],NROWS(TestMatrix_1_2x3),NCOL(TestMatrix_1_2x3));

    mtx_init_f64(&UkfMat.Wm_weight_vector, &Wm_sigma_weight_1x9[0][0],NROWS(Wm_sigma_weight_1x9),NCOL(Wm_sigma_weight_1x9));

    mtx_init_f64(&UkfMat.Wc_weight_vector, &Wc_sigma_weight_1x9[0][0],NROWS(Wc_sigma_weight_1x9),NCOL(Wc_sigma_weight_1x9));

    mtx_init_f64(&UkfMat.u_system_input, &u_curr_system_input_4x1[0][0],NROWS(u_curr_system_input_4x1),NCOL(u_curr_system_input_4x1));

    mtx_init_f64(&UkfMat.u_prev_system_input, &u_prev_system_input_4x1[0][0],NROWS(u_prev_system_input_4x1),NCOL(u_prev_system_input_4x1));

    mtx_init_f64(&UkfMat.y_meas, &y_curr_system_meas_2x1[0][0],NROWS(y_curr_system_meas_2x1),NCOL(y_curr_system_meas_2x1));

    mtx_init_f64(&UkfMat.y_predicted_mean, &y_mean_system_predict_2x1[0][0],NROWS(y_mean_system_predict_2x1),NCOL(y_mean_system_predict_2x1));

    mtx_init_f64(&UkfMat.x_system_states, &x_system_states_4x1[0][0],NROWS(x_system_states_4x1),NCOL(x_system_states_4x1));

    mtx_init_f64(&UkfMat.x_system_states_ic, &x_system_states_ic_4x1[0][0],NROWS(x_system_states_ic_4x1),NCOL(x_system_states_ic_4x1));

    mtx_init_f64(&UkfMat.x_system_states_correction, &x_system_states_correction_4x1[0][0],NROWS(x_system_states_correction_4x1),NCOL(x_system_states_correction_4x1));

    mtx_init_f64(&UkfMat.X_sigma_points, &X_sigma_points_4x9[0][0],NROWS(X_sigma_points_4x9),NCOL(X_sigma_points_4x9));

    mtx_init_f64(&UkfMat.Y_sigma_points, &Y_sigma_points_2x9[0][0],NROWS(Y_sigma_points_2x9),NCOL(Y_sigma_points_2x9));

    mtx_init_f64(&UkfMat.Pxx_error_covariance, &Px_state_cov_4x4[0][0],NROWS(Px_state_cov_4x4),NCOL(Px_state_cov_4x4));

    mtx_init_f64(&UkfMat.Pxx_covariance_correction, &Pxx_covariance_correction_4x4[0][0],NROWS(Pxx_covariance_correction_4x4),NCOL(Pxx_covariance_correction_4x4));

    mtx_init_f64(&UkfMat.Pyy_out_covariance, &Pyy_out_cov_2x2[0][0],NROWS(Pyy_out_cov_2x2),NCOL(Pyy_out_cov_2x2));

    mtx_init_f64(&UkfMat.Pyy_out_covariance_copy, &Pyy_out_cov_copy_2x2[0][0],NROWS(Pyy_out_cov_copy_2x2),NCOL(Pyy_out_cov_copy_2x2));

    mtx_init_f64(&UkfMat.Ryy0_init_out_covariance, &Ryy_out_cov_noise_2x2[0][0],NROWS(Ryy_out_cov_noise_2x2),NCOL(Ryy_out_cov_noise_2x2));

    mtx_init_f64(&UkfMat.Pxy_cross_covariance, &Pxy_state_out_cov_4x2[0][0],NROWS(Pxy_state_out_cov_4x2),NCOL(Pxy_state_out_cov_4x2));

    mtx_init_f64(&UkfMat.K_kalman_gain, &K_kalman_gain_4x2[0][0],NROWS(K_kalman_gain_4x2),NCOL(K_kalman_gain_4x2));

    mtx_init_f64(&UkfMat.K_kalman_gain_transp, &K_kalman_transp_gain_2x4[0][0],NROWS(K_kalman_transp_gain_2x4),NCOL(K_kalman_transp_gain_2x4));

    mtx_init_f64(&UkfMat.I_identity_matrix, &temporal_2x2[0][0],NROWS(temporal_2x2),NCOL(temporal_2x2));

    mtx_init_f64(&UkfMat.Qxx_process_noise_cov, &Qxx_process_noise_cov_4x4[0][0],NROWS(Qxx_process_noise_cov_4x4),NCOL(Qxx_process_noise_cov_4x4));

    mtx_init_f64(&UkfMat.Pxx0_init_error_covariance, &P0_state_cov_4x4[0][0],NROWS(P0_state_cov_4x4),NCOL(P0_state_cov_4x4));//

    UkfMat.fcnPredict = &PredictFcn[0];
    UkfMat.fcnObserve = &ObservFcn[0];

    tfInitFail = ukf_init(&ukfIo,param,4,2, &UkfMat);

    if(tfInitFail == 0)
    {             
        for(simLoop=1;simLoop<15;simLoop++)
        {
            
            y_curr_system_meas_2x1[0][0] = yt[0][simLoop];
            y_curr_system_meas_2x1[1][0] = yt[1][simLoop];
            
            (void)ukf_step(&ukfIo);
            
            //printf("system states \n");
            printf("%2.14f  %2.14f  %2.14f  %2.14f ", x_system_states_4x1[0][0], x_system_states_4x1[1][0], x_system_states_4x1[2][0], x_system_states_4x1[3][0]);
            //show_matrix_obj(UkfMat.x_system_states);
            printf("\n");
            printf("%2.14f  %2.14f  %2.14f  %2.14f ", x_system_states_4x1[0][0]-x_exp[simLoop-1][0], x_system_states_4x1[1][0]-x_exp[simLoop-1][1],x_system_states_4x1[2][0]-x_exp[simLoop-1][2], x_system_states_4x1[3][0]-x_exp[simLoop-1][3]);
            printf("\n");
            
            rootSquareErr_X0 += fabs(x_system_states_4x1[0][0]-x_exp[simLoop-1][0]);
            rootSquareErr_X1 += fabs(x_system_states_4x1[1][0]-x_exp[simLoop-1][1]);
            rootSquareErr_X2 += fabs(x_system_states_4x1[2][0]-x_exp[simLoop-1][2]);
            rootSquareErr_X3 += fabs(x_system_states_4x1[3][0]-x_exp[simLoop-1][3]);
        }
        printf("\n");
        printf("%2.16f  %2.16f  %2.16f  %2.16f ",rootSquareErr_X0, rootSquareErr_X1,rootSquareErr_X2,rootSquareErr_X3);
    }
    else
    {
        //initialization fail
        //TBD
    }
}
//




