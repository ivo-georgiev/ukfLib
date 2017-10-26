#ifndef MTXLIB_FILE
#define MTXLIB_FILE

#include "System_Types.h"
#include<math.h>
/*---------------------------------------------*/
/*         Macros definiton                    */
/*---------------------------------------------*/
#define NCOL(arr) sizeof(arr[0])/sizeof(arr[0][0])
#define NROWS(arr)    sizeof(arr)/sizeof(arr[0])
//#define rowxcol(arr) sizeof(arr)/sizeof(arr[0][0])


int mtx_chol1_f64(double* mtxA, double* mtxL,uint8 mtxSize);

typedef int mtxResultInfo;
#define MTX_OPERATION_OK (mtxResultInfo)0

#define MTX_SINGULAR (mtxResultInfo)251
#define MTX_SIZE_MISMATCH (mtxResultInfo)252
#define MTX_NOT_SQUARE (mtxResultInfo)253
#define MTX_NOT_POS_DEFINED (mtxResultInfo)254
#define MTX_OPERATION_ERROR (mtxResultInfo)255

typedef struct tMatrix
{
    uint8 nrow;
    uint8 ncol;
    double* val;
}tMatrix;

extern mtxResultInfo mtx_init_f64(tMatrix* A, double * M, uint8 nrow, uint8 ncol);
extern mtxResultInfo mtx_mul_f64(tMatrix * pA, tMatrix * pB, tMatrix * pC);
extern mtxResultInfo mtx_transp_square_f64(tMatrix * pA);
extern mtxResultInfo mtx_transp_dest_f64(tMatrix * pA,tMatrix * pB);
extern mtxResultInfo mtx_diagsum_f64(tMatrix * pA, double * diagsum);
extern mtxResultInfo mtx_chol_upper_f64(tMatrix * pA);
extern mtxResultInfo mtx_chol_lower_f64(tMatrix * pA);
extern mtxResultInfo mtx_inv_f64(tMatrix * pA, tMatrix * pI);
extern mtxResultInfo mtx_add_f64(tMatrix * pA,tMatrix * pB);
extern mtxResultInfo mtx_subtract_f64(tMatrix * pA,tMatrix * pB);
extern mtxResultInfo mtx_mul_scalar_f64(tMatrix * pA,double scalar);
extern mtxResultInfo mtx_add_scalar_f64(tMatrix * pA,double scalar);
extern mtxResultInfo mtx_subtract_scalar_f64(tMatrix * pA,double scalar);
extern mtxResultInfo mtx_cpy_f64(tMatrix * const pDestP,tMatrix const * const pSrcP);
extern mtxResultInfo mtx_identity_f64(tMatrix * pI);
extern mtxResultInfo mtx_zeros_f64(tMatrix * pZ);
//old integerlib
#endif

