#include "../cfg/ukfCfg.h"
#include "../cfg/ukfCfg1.h"

#define cfg0 (uint8)0
#define cfg1 (uint8)1

/*---------------------------------------------*/
/*           Function Prototype                */
/*---------------------------------------------*/
void show_matrix_obj(tMatrix A);
void show_matrix(float64 * A, int n,int m);
void ukf_test(void);
void mtxlib_test(void);
void report_compiler(void);

int main(void)
{
	report_compiler();

    //generic matrix operation test
    mtxlib_test();

    //UKF test start here
    ukf_test();

    return 0;
}

void show_matrix_obj(tMatrix A)
{
    int i,j;
    
    for (i = 0; i < A.nrow; i++) 
    {
        for (j = 0; j < A.ncol; j++)
        {
            printf("%2.14f ", A.val[A.ncol*i + j]);
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
            printf("%2.14f ", A[m*i + j]);
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
    _Bool tfInitCfg0 = 0;
    _Bool tfInitCfg1 = 0;
    tUKF ukfIo[2];
    uint32 simLoop;

    //UKF filter measurement input(data log is generated in matlab and used for UKF simulation for 15 iteration) 
    static const float64 yt[2][15]=
    {
        {0,  16.085992708563385,  12.714829185978214,  14.528500994457660,  19.105561355310275, 23.252820029388918,  29.282949862903255,  36.270058819651275,  44.244884173240955,  47.394243121124411, 55.988905459180458, 61.667450941562109, 68.624980301613647, 76.337963872393104, 82.611325690835159}, //y1 test
        {0,  16.750821420874981,  14.277640835870006,  16.320754051600520,  20.560460303503849, 24.827446289454556,  31.290961393448615,  36.853553457560210,  42.157283183453522,  49.382835230961490, 57.516319669684677, 65.664496283509095, 71.428712755732704, 79.241720894223079, 84.902760328915676 }
    };
    
    //UKF filter expected system states calculated with matlab script for 15 iterations
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

    //UKF initialization: CFG0
    tfInitCfg0 = ukf_init(&ukfIo[cfg0], &UkfMatrixCfg0);

    if(tfInitCfg0 == 0 )
    {   
        float64 err[4]={0,0,0,0};
        float64 absErrAccum[4] = {0,0,0,0};

        //UKF simulation CFG0: BEGIN
        for(simLoop=1;simLoop<15;simLoop++)
        {
            float64 * const py_cfg0 = ukfIo[cfg0].input.y.val;

            //UKF:CFG0 apply/load system measurements in working array for current iteration.
            py_cfg0[0] = yt[0][simLoop];
            py_cfg0[1] = yt[1][simLoop];
            
            //UKF:CFG0 periodic task call
            (void)ukf_step(&ukfIo[cfg0]);
            
            err[0] = fabs(ukfIo[cfg0].update.x.val[0] - x_exp[simLoop-1][0]);
            err[1] = fabs(ukfIo[cfg0].update.x.val[1] - x_exp[simLoop-1][1]);
            err[2] = fabs(ukfIo[cfg0].update.x.val[2] - x_exp[simLoop-1][2]);
            err[3] = fabs(ukfIo[cfg0].update.x.val[3] - x_exp[simLoop-1][3]);
            
            printf("Loop: %d |system states : ukf.m | system states : est | system states : impl. diff \n",(int)simLoop);
            printf("          %2.14f        %2.14f       %2.14f\n", x_exp[simLoop-1][0], ukfIo[cfg0].update.x.val[0], err[0]);
            printf("          %2.14f        %2.14f       %2.14f\n", x_exp[simLoop-1][1], ukfIo[cfg0].update.x.val[1], err[1]);
            printf("          %2.14f        %2.14f       %2.14f\n", x_exp[simLoop-1][2], ukfIo[cfg0].update.x.val[2], err[2]);
            printf("          %2.14f        %2.14f       %2.14f\n", x_exp[simLoop-1][3], ukfIo[cfg0].update.x.val[3], err[3]);           
            
            //accumulate the differennce between reference matlab implementation and results from C code execution 
            absErrAccum[0] += err[0];
            absErrAccum[1] += err[1];
            absErrAccum[2] += err[2];
            absErrAccum[3] += err[3];
        }
        printf("Accumulated error: CFG0 \n");
        printf("%2.16f  \n%2.16f  \n%2.16f  \n%2.16f \n",absErrAccum[0], absErrAccum[1],absErrAccum[2],absErrAccum[3]);

        //UKF simulation CFG0: END
    }
    else
    {
        //initialization fail
    }

    //UKF initialization: CFG1(free pendulum)
    tfInitCfg1 = ukf_init(&ukfIo[cfg1], &UkfMatrixCfg1);

    if(tfInitCfg1 == 0 )
    {        
        static float64 tetha = 0.5;  //initial conditions for angle
        static float64 tetha_dot = 0;//initial conditions for angle speed
        const float64 B = 0.05; //kg*s/m 
        const float64 l = 0.613;
        const float64 m = 0.5;
        const float64 g = 9.81;
        static const float64 T0 = 0.0001;

        //UKF simulation: BEGIN
        for(simLoop=0;simLoop<70;simLoop++)
        {
            float64 * const py_cfg1 = ukfIo[cfg1].input.y.val;
            float64 err[2]={0,0};
            
            //UKF:CFG1 apply/load system measurements in working array for current iteration
                 
            tetha = tetha + T0*tetha_dot;
            tetha_dot = tetha_dot - ((T0*B*tetha_dot)/m) - ((T0*g)/l)*sin(tetha);

            py_cfg1[0] = tetha;
            
            //UKF:CFG0 periodic task call
            (void)ukf_step(&ukfIo[cfg1]);			  
            
            err[0] = fabs(ukfIo[cfg1].update.x.val[0] - tetha);
            err[1] = fabs(ukfIo[cfg1].update.x.val[1] - tetha_dot);         
            
            printf("Loop: %d |system states : real | system states : est | system states : err \n",(int)simLoop);
            printf("          %2.14f       %2.14f      %2.14f\n", tetha, ukfIo[1].update.x.val[0], err[0]);
            printf("          %2.14f      %2.14f     %2.14f\n", tetha_dot, ukfIo[1].update.x.val[1], err[1]);
            
        }           
        //UKF simulation: END
    }
    else
    {
        //initialization fail
        //TBD
    }
}
/******************************************************************************************************************************************************************************************************\
 ***  FUNCTION:
 ***      void mtxlib_test(void)
 *** 
 ***  DESCRIPTION:
 ***      Test some generic matrix operations from mtxLib.c      
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
void mtxlib_test(void)
{
    /*------------------------------------------------*/
    /*             Sample Matrix not related with UKF */
    /*------------------------------------------------*/
    static float64 symMtx[5][5] =
    {{ 3.009964742796415,  -1.009481719676593,   1.161774056277037,  -0.291928007896218,  -0.775215823251770},
    { -1.009481719676593,   3.009964742796415,  -1.009481719676593,   1.161774056277037,  -0.291928007896218},
    {  1.161774056277037,  -1.009481719676593,   3.009964742796415,  -1.009481719676593,   1.161774056277037},
    { -0.291928007896218,   1.161774056277037,  -1.009481719676593,   3.009964742796415,  -1.009481719676593},
    { -0.775215823251770,  -0.291928007896218,   1.161774056277037,  -1.009481719676593,   3.009964742796415}};
    
    static float64 symMtxChol[5][5] =
    {{1.734924996302842,     -0.581858997840148,   0.669639355449256,      -0.168265491890613,  -0.446829589119858},
    {0                 ,      1.634443284249678,  -0.379239855780691,       0.650904975441148,  -0.337680621986338},
    {0                 ,      0                ,   1.554903536627710,      -0.418003689501540,   0.857240820764834},
    {0                 ,      0                ,   0                ,       1.543776893059448,  -0.328117294491480},
    {0                 ,      0                ,   0                ,       0                ,   1.361527478565284}};
    
    static float64 Identity_5x5[4][4] =
    {{1.0,  0     ,0,     0,     },
    {0,     1.0   ,0,     0,     },
    {0,     0     ,1.0,   0,     },
    {0,     0     ,0,     1.0   }};
    
    static float64 TestMatrix_0_4x4[4][4] =
    {{3.0,  5.0, -1.0,  -4},
    { 1.0,  4.0, -0.7,  -3},
    { 0,   -2.0,  0,     1},
    {-2.0,  6.0,  0,     0.3}};
    
    
#if 0
    static float64 TestMatrix_1_3x3[3][3] = 
    {{10.5, 2.17, 3.03},
    { 0.44, 0.59, 6.89},
    { 7.56, 8.17, 9.21}};
    
    static float64 TestMatrix_2_3x3[3][3] = 
    { {1.11, 29.3, 31.2},
    {45.3, 5.17, 6.11},
    {7.61, 88.0, 9.34}};
#endif
    
    static float64 TestMatrix_1_2x3[2][3] = 
    { {1.11, 29.3, 31.2},   //size 3x3
    {45.3, 5.17, 6.11}};
    
    static float64 TestMatrixDest_3x2[3][2] =
    {
        {0,0},
        {0,0},
        {0,0}  
    };
    tMatrix myFactMatrix;
    tMatrix myTestMatx={0,0,0,NULL};
    tMatrix myChol={0,0,0,NULL};
    tMatrix Im={0,0,0,NULL};
    tMatrix oTestMatrix_0_4x4={0,0,0,NULL};
    tMatrix oTestMatrixDest_3x2={0,0,0,NULL};
    
    mtx_init_f64(&myTestMatx,&TestMatrix_1_2x3[0][0],NROWS(TestMatrix_1_2x3),NCOL(TestMatrix_1_2x3),COLXROW(TestMatrix_1_2x3));
    mtx_init_f64(&myChol,&symMtxChol[0][0],NROWS(symMtxChol),NCOL(symMtxChol),COLXROW(symMtxChol));
    mtx_init_f64(&Im,&Identity_5x5[0][0],NROWS(Identity_5x5),NCOL(Identity_5x5),COLXROW(Identity_5x5));
    mtx_init_f64(&oTestMatrix_0_4x4,&TestMatrix_0_4x4[0][0],NROWS(TestMatrix_0_4x4),NCOL(TestMatrix_0_4x4),COLXROW(TestMatrix_0_4x4));
    mtx_init_f64(&oTestMatrixDest_3x2,&TestMatrixDest_3x2[0][0],NROWS(TestMatrixDest_3x2),NCOL(TestMatrixDest_3x2),COLXROW(TestMatrixDest_3x2));
    //show_matrix(&TestMatrix_1_2x3[0][0],2,3);
    //show_matrix_obj(myTestMatx);
    
    (void)mtx_init_f64(&myFactMatrix,&symMtx[0][0],NROWS(symMtx),NCOL(symMtx),COLXROW(symMtx));
    //show_matrix(&symMtx[0][0],5,5);
    
    (void)mtx_chol_lower_f64(&myFactMatrix);
    show_matrix_obj(myFactMatrix);
    
    /*show_matrix_obj(myChol);
    mtx_transp_square_f64(&myChol);
    show_matrix_obj(myChol);
    */
    
    show_matrix_obj(Im);
    show_matrix_obj(oTestMatrix_0_4x4);
    mtx_inv_f64(&oTestMatrix_0_4x4, &Im);
    show_matrix_obj(Im);

    mtx_identity_f64(&Im);
    show_matrix_obj(Im);
   
    //show_matrix_obj(myTestMatx);
    //mtx_transp_dest_f64(&myTestMatx,&oTestMatrixDest_3x2);
    //show_matrix_obj(oTestMatrixDest_3x2);*/
}

void report_compiler(void)
{
	fprintf(stderr, "sizeof float = %d\nsizeof double = %d\nsizeof long double = %d\n", 8 * __SIZEOF_FLOAT__, 8 *  __SIZEOF_DOUBLE__, 8 * __SIZEOF_LONG_DOUBLE__);
}

