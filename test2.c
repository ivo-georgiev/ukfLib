#include<stdio.h>
#include<math.h>
#include <stdlib.h>
#include "mtxLib.h"
#include "ukfLib.h"
#include "ukfCfg.h"

/*---------------------------------------------*/
/*           Function Prototype                */
/*---------------------------------------------*/
void show_matrix_obj(tMatrix A);
void show_matrix(double * A, int n,int m);
void ukf_test(void);

//


/*---------------------------------------------*/
/*             Sample Matrix                   */
/*---------------------------------------------*/
double symMtx[5][5] =
{{ 3.009964742796415,  -1.009481719676593,   1.161774056277037,  -0.291928007896218,  -0.775215823251770},
{ -1.009481719676593,   3.009964742796415,  -1.009481719676593,   1.161774056277037,  -0.291928007896218},
{  1.161774056277037,  -1.009481719676593,   3.009964742796415,  -1.009481719676593,   1.161774056277037},
{ -0.291928007896218,   1.161774056277037,  -1.009481719676593,   3.009964742796415,  -1.009481719676593},
{ -0.775215823251770,  -0.291928007896218,   1.161774056277037,  -1.009481719676593,   3.009964742796415}};

double symMtxChol[5][5] =
{{1.734924996302842,     -0.581858997840148,   0.669639355449256,      -0.168265491890613,  -0.446829589119858},
{0                 ,      1.634443284249678,  -0.379239855780691,       0.650904975441148,  -0.337680621986338},
{0                 ,      0                ,   1.554903536627710,      -0.418003689501540,   0.857240820764834},
{0                 ,      0                ,   0                ,       1.543776893059448,  -0.328117294491480},
{0                 ,      0                ,   0                ,       0                ,   1.361527478565284}};

double Identity_5x5[4][4] =
{{1.0,  0     ,0,     0,     },
{0,     1.0   ,0,     0,     },
{0,     0     ,1.0,   0,     },
{0,     0     ,0,     1.0   }};

double TestMatrix_0_4x4[4][4] =
{{3.0,  5.0, -1.0,  -4},
{ 1.0,  4.0, -0.7,  -3},
{ 0,   -2.0,  0,     1},
{-2.0,  6.0,  0,     0.3}};


double TestMatrix_1_3x3[3][3] = 
{{10.5, 2.17, 3.03},
{ 0.44, 0.59, 6.89},
{ 7.56, 8.17, 9.21}};

double TestMatrix_2_3x3[3][3] = 
{ {1.11, 29.3, 31.2},
{45.3, 5.17, 6.11},
{7.61, 88.0, 9.34}};


double TestMatrix_1_2x3[2][3] = 
{ {1.11, 29.3, 31.2},   //size 3x3
  {45.3, 5.17, 6.11}};

double TestMatrixDest_3x2[3][2] =
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

void show_matrix(double * A, int n,int m)
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

int MtxSumDiag(int*A ,int m, int n)
{
    int i,j,sum=0;
    for(i=0;i<m;i++)
    {
        for(j=0;j<n;j++)
        {
            if(i == j)
            {
                sum += *(A+i*n+j);
            }
        }
    } 
    return sum;
}



void ukf_test(void)
{
    double param[3]={1,2,0};
    tUkfMatrix UkfMat;
    tUKF ukfIo;
    tMatrix myTestMatx;
    int simLoop;

    static const double yt[2][5]=
    {
        {15.025369860836, 12.705439388362, 15.965607350939, 19.788614274432, 27.715926478831}, //y1 test
        {15.413169415429, 14.252524720158, 15.482947630858, 19.986301110182, 22.473922714307}
    };

    static const double xt[4][5]=
    {
        {0.350098526093, 5.368706332304, 9.935544766314, 14.415664795863, 18.625807437381 },
        {-0.291307157467, 4.766352787570, 9.996447396680, 15.295066899099, 21.200671094783 },
        {50.186078062109, 45.668384340102, 44.801200295491, 42.101426415178, 43.530912222830},
        {50.576599450369, 52.300946091105, 52.986195024183, 59.056041956846, 58.646109840247 }

    };



    mtx_init_f64(&myTestMatx,&TestMatrix_1_2x3[0][0],NROWS(TestMatrix_1_2x3),NCOL(TestMatrix_1_2x3));

    mtx_init_f64(&UkfMat.Wm_weight_vector, &Wm_sigma_weight_1x9[0],1,9);

    mtx_init_f64(&UkfMat.Wc_weight_vector, &Wc_sigma_weight_1x9[0],1,9);

    mtx_init_f64(&UkfMat.u_system_input, &u_curr_system_input_4x1[0],1,4);

    mtx_init_f64(&UkfMat.u_prev_system_input, &u_prev_system_input_4x1[0],1,4);

    mtx_init_f64(&UkfMat.y_meas, &y_curr_system_meas_2x1[0],2,1);

    mtx_init_f64(&UkfMat.y_predicted_mean, &y_mean_system_predict_2x1[0],2,1);

    mtx_init_f64(&UkfMat.x_system_states, &x_system_states_4x1[0],4,1);

    mtx_init_f64(&UkfMat.x_system_states_correction, &x_system_states_correction_4x1[0],4,1);

    mtx_init_f64(&UkfMat.X_sigma_points, &X_sigma_points_4x9[0][0],4,9);

    mtx_init_f64(&UkfMat.Y_sigma_points, &Y_sigma_points_2x9[0][0],2,9);

    mtx_init_f64(&UkfMat.Pxx_error_covariance, &Px_state_cov_4x4[0][0],4,4);

    mtx_init_f64(&UkfMat.Pxx_covariance_correction, &Pxx_covariance_correction_4x4[0][0],4,4);

    mtx_init_f64(&UkfMat.Pyy_out_covariance, &Pyy_out_cov_2x2[0][0],2,2);

    mtx_init_f64(&UkfMat.Ryy0_init_out_covariance, &Ryy_out_cov_noise_2x2[0][0],2,2);

    mtx_init_f64(&UkfMat.Pxy_cross_covariance, &Pxy_state_out_cov_4x2[0][0],4,2);

    mtx_init_f64(&UkfMat.K_kalman_gain, &K_kalman_gain_4x2[0][0],4,2);

    mtx_init_f64(&UkfMat.K_kalman_gain_transp, &K_kalman_transp_gain_2x4[0][0],2,4);

    mtx_init_f64(&UkfMat.I_identity_matrix, &temporal_2x2[0][0],2,2);

    mtx_init_f64(&UkfMat.Qxx_process_noise_cov, &Qxx_process_noise_cov_4x4[0][0],4,4);

    mtx_init_f64(&UkfMat.Pxx0_init_error_covariance, &P0_state_cov_4x4[0][0],4,4);

    UkfMat.fcnPredict = &PredictFcn[0];
    UkfMat.fcnObserve = &ObservFcn[0];

    (void)ukf_init(&ukfIo,param,4,2, &UkfMat);

    for(simLoop=0;simLoop<5;simLoop++)
    {
        y_curr_system_meas_2x1[0] = yt[0][simLoop];
        y_curr_system_meas_2x1[1] = yt[1][simLoop];

        (void)ukf_step(&ukfIo);
    }
   

}
//




