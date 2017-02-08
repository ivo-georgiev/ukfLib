


int mtxLib_Cholesky1_LL_dp(double* mtxA, double* mtxL,int mtxSize);

typedef int mtxResultInfo;
#define MTX_OPERATION_OK (mtxResultInfo)0

#define MTX_NOT_POS_DEFINED (mtxResultInfo)254
#define MTX_OPERATION_ERROR (mtxResultInfo)255

typedef struct sMatrixType
{
    int nrow;
    int ncol;
    double* val;
}sMatrixType;

mtxResultInfo mtxLib_Mul_dp(sMatrixType A, sMatrixType B, sMatrixType C);
mtxResultInfo mtxLib_Transp_dp(sMatrixType A);
mtxResultInfo mtx_diagsum_dp(sMatrixType A, double * DiagSumP);
mtxResultInfo mtxLib_Cholesky_LL_dp(sMatrixType A);

//old integerlib
