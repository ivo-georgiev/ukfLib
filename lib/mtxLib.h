/******************************************************************************************************************************************************************************************************\
 ***
 *** Description       : IMPLEMENTATION OF BASIC MATRIX OPERATION
 *** Codefile          : mtxLib.h
 ***
 *** MIT License
 ***
 *** Copyright (c) 2017 ivo-georgiev
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
#ifndef MTXLIB_FILE
#define MTXLIB_FILE

#include "System_Types.h"
#include <math.h>
#include <stdlib.h>
/*---------------------------------------------*/
/*         Macros definiton                    */
/*---------------------------------------------*/
#define NCOL(arr) sizeof(arr[0]) / sizeof(arr[0][0])
#define NROWS(arr) sizeof(arr) / sizeof(arr[0])
#define COLXROW(arr) sizeof(arr) / sizeof(arr[0][0])

enum mtxResultInfo mtx_chol1(double *A, double *L, ptrdiff_t size);

enum mtxResultInfo
{
	MTX_OPERATION_OK = 0,
	MTX_SINGULAR = 251,
	MTX_SIZE_MISMATCH = 252,
	MTX_NOT_SQUARE = 253,
	MTX_NOT_POS_DEFINED = 254,
	MTX_OPERATION_ERROR = 255,
};

typedef struct sMatrix
{
	ptrdiff_t nelem;
	ptrdiff_t nrow;
	ptrdiff_t ncol;
	double *val;
} Matrix_t;

typedef struct sMatrixBool
{
	ptrdiff_t nelem;
	ptrdiff_t nrow;
	ptrdiff_t ncol;
	_Bool *val;
} MatrixBool_t;

extern enum mtxResultInfo mtx_init_bool(MatrixBool_t *const pSrc, _Bool *const pValue, const ptrdiff_t nrow, const ptrdiff_t ncol, const ptrdiff_t nelem);
extern enum mtxResultInfo mtx_init(Matrix_t *const pSrc, double *const pValue, const ptrdiff_t nrow, const ptrdiff_t ncol, const ptrdiff_t nelem);
extern enum mtxResultInfo mtx_mul(Matrix_t const *const pSrc1, Matrix_t const *const pSrc2, Matrix_t *const pDst);
extern enum mtxResultInfo mtx_transp_square(Matrix_t *const pSrc);
extern enum mtxResultInfo mtx_transp_dest(Matrix_t const *const pSrc, Matrix_t *const pDst);
extern enum mtxResultInfo mtx_diagsum(Matrix_t *pSrc, double *diagsum);
extern enum mtxResultInfo mtx_chol_upper(Matrix_t *const pSrc);
extern enum mtxResultInfo mtx_chol_lower(Matrix_t *const pSrc);
extern enum mtxResultInfo mtx_inv(Matrix_t *const pSrc, Matrix_t *const pDst);
extern enum mtxResultInfo mtx_add(Matrix_t *const pDst, Matrix_t const *const pSrc);
extern enum mtxResultInfo mtx_sub(Matrix_t *const pDst, Matrix_t const *const pSrc);
extern enum mtxResultInfo mtx_mul_scalar(Matrix_t *const pSrc, const double scalar);
extern enum mtxResultInfo mtx_add_scalar(Matrix_t *const pSrc, const double scalar);
extern enum mtxResultInfo mtx_sub_scalar(Matrix_t *const pSrc, const double scalar);
extern enum mtxResultInfo mtx_cpy(Matrix_t *const pDst, Matrix_t const *const pSrc);
extern enum mtxResultInfo mtx_identity(Matrix_t *const pSrc);
extern enum mtxResultInfo mtx_zeros(Matrix_t *const pSrc);
extern enum mtxResultInfo mtx_mul_src2tr(Matrix_t const *const pSrc1, Matrix_t const *const pSrc2, Matrix_t *const pDst);

#endif
