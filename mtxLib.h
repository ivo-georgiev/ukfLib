/*---------------------------------------------*/
/*         Macros definiton                    */
/*---------------------------------------------*/
#define NCOL(arr) sizeof(arr[0])/sizeof(arr[0][0])
#define NROWS(arr)    sizeof(arr)/sizeof(arr[0])
//#define rowxcol(arr) sizeof(arr)/sizeof(arr[0][0])


int mtx_chol1_dp(double* mtxA, double* mtxL,int mtxSize);

typedef int mtxResultInfo;
#define MTX_OPERATION_OK (mtxResultInfo)0

#define MTX_DIMENSION_ERR (mtxResultInfo)252
#define MTX_NOT_SQUARE (mtxResultInfo)253
#define MTX_NOT_POS_DEFINED (mtxResultInfo)254
#define MTX_OPERATION_ERROR (mtxResultInfo)255

typedef struct sMatrixType
{
    int nrow;
    int ncol;
    double* val;
}sMatrixType;

mtxResultInfo mtx_init_dp(sMatrixType* A, double * M, int nrow, int ncol);
mtxResultInfo mtx_mul_dp(sMatrixType A, sMatrixType B, sMatrixType C);
mtxResultInfo mtx_transp_dp(sMatrixType A);
mtxResultInfo mtx_diagsum_dp(sMatrixType A, double * diagsum);
mtxResultInfo mtx_chol_dp(sMatrixType A);

//old integerlib
