/******************************************************************************************************************************************************************************************************\
 ***
 *** Description       : IMPLEMENTATION OF BASIC MATRIX OPERATION
 *** Codefile          : mtxLib.c
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

#include "mtxLib.h"

enum mtxResultInfo mtx_init_bool(MatrixBool_t *const pSrc, _Bool *const pValue, const uint8_t nrow, const uint8_t ncol, const uint16_t nelem)
{
	pSrc->val = pValue;
	pSrc->ncol = ncol;
	pSrc->nrow = nrow;
	pSrc->nelem = nelem;
	return MTX_OPERATION_OK;
}
enum mtxResultInfo mtx_init(Matrix_t *const pSrc, double *const pValue, const uint8_t nrow, const uint8_t ncol, const uint16_t nelem)
{
	pSrc->val = pValue;
	pSrc->ncol = ncol;
	pSrc->nrow = nrow;
	pSrc->nelem = nelem;
	return MTX_OPERATION_OK;
}
/**
 * For square matrix only
 */
enum mtxResultInfo mtx_diagsum(Matrix_t *pSrc, double *diagsum)
{
	enum mtxResultInfo Result = MTX_OPERATION_OK;
	double const *const pSrcL = (double *)pSrc->val;
	const uint8_t ncol = pSrc->ncol;
	uint16_t eIdx;
	double sum = pSrcL[0];

	if (pSrc->nrow == ncol)
	{
		for (eIdx = 1; eIdx < pSrc->nelem; eIdx++)
		{
			const uint16_t cmpLeft = (uint16_t)(eIdx / ncol);

			sum += eIdx < ncol ? 0 : cmpLeft == eIdx % (cmpLeft * ncol) ? pSrcL[eIdx] : 0;
		}
	}
	else
	{
		Result = MTX_SIZE_MISMATCH;
	}

	*diagsum = sum;

	return Result;
}
/**
 * @brief  \f$A=A'\f$ or \f$A=A^{\tau}\f$
 */
enum mtxResultInfo mtx_transp_square(Matrix_t *const pSrc)
{
	enum mtxResultInfo ResultL = MTX_OPERATION_OK;
	const uint8_t nrow = pSrc->nrow;
	const uint8_t ncol = pSrc->ncol;
	double *const pSrcL = (double *)pSrc->val;
	uint8_t row, col;
	double temp;

	if (nrow == ncol)
	{
		for (row = 0; row < nrow; row++)
		{
			for (col = 0; col < ncol; col++)
			{
				if (row != col && row < col)
				{
					temp = pSrcL[nrow * row + col];
					pSrcL[ncol * row + col] = pSrcL[ncol * col + row];
					pSrcL[ncol * col + row] = temp;
				}
			}
		}
	}
	else
	{
		ResultL = MTX_NOT_SQUARE;
	}

	return ResultL;
}
enum mtxResultInfo mtx_transp_dest(Matrix_t const *const pSrc, Matrix_t *const pDst)
{
	enum mtxResultInfo ResultL = MTX_OPERATION_OK;
	double const *const pSrcL = (double *)pSrc->val;
	double *const pDstL = (double *)pDst->val;
	const uint8_t nRowSrcL = pSrc->nrow;
	const uint8_t nColSrcL = pSrc->ncol;
	const uint8_t nRowDstL = pDst->nrow;
	const uint8_t nColDstL = pDst->ncol;
	uint8_t row, col;

	if (nRowSrcL == nColDstL || nColSrcL == nRowDstL)
	{
		for (row = 0; row < nRowDstL; row++)
		{
			for (col = 0; col < nColDstL; col++)
			{
				pDstL[nColDstL * row + col] = pSrcL[nColSrcL * col + row];
			}
		}
	}
	else
	{
		ResultL = MTX_SIZE_MISMATCH;
	}

	return ResultL;
}
/**
 * @brief \f$C=A \cdot B\f$
 *
 * Matrix multiplication
 *
 */
enum mtxResultInfo mtx_mul(Matrix_t const *const pSrc1, Matrix_t const *const pSrc2, Matrix_t *const pDst)
{
	enum mtxResultInfo ResultL = MTX_OPERATION_OK;
	double const *const pSrc1L = (double *)pSrc1->val;
	double const *const pSrc2L = (double *)pSrc2->val;
	double *const pDstL = (double *)pDst->val;
	uint8_t row, col, k;
	double sum;

	if (pSrc1->ncol == pSrc2->nrow)
	{
		for (row = 0; row < pSrc1->nrow; row++)
		{
			for (col = 0; col < pSrc2->ncol; col++)
			{
				sum = 0;
				for (k = 0; k < pSrc1->ncol; k++)
				{
					sum += pSrc1L[pSrc1->ncol * row + k] * pSrc2L[pSrc2->ncol * k + col];
				}
				pDstL[pDst->ncol * row + col] = sum;
			}
		}
	}
	else
	{
		ResultL = MTX_SIZE_MISMATCH;
	}

	return ResultL;
}
/**
 * @brief Special multiplication \f$Dst=Src1 \cdot Src2^{\tau}\f$
 *
 *  Special function for multiplication of Src1 matrix with transpose image of
 *  Src2 matrix. Be sure that array size are suitable for multiplication after
 *  Src2 transpose
 */
enum mtxResultInfo mtx_mul_src2tr(Matrix_t const *const pSrc1, Matrix_t const *const pSrc2, Matrix_t *const pDst)
{
	enum mtxResultInfo ResultL = MTX_OPERATION_OK;
	double const *const pSrc1L = (double *)pSrc1->val;
	double const *const pSrc2L = (double *)pSrc2->val;
	double *const pDstL = (double *)pDst->val;
	uint8_t rowSrc1, rowSrc2, k;
	double sum;

	if (pSrc1->ncol == pSrc2->ncol)
	{
		for (rowSrc1 = 0; rowSrc1 < pSrc1->nrow; rowSrc1++)
		{
			for (rowSrc2 = 0; rowSrc2 < pSrc2->nrow; rowSrc2++)
			{
				sum = 0;
				for (k = 0; k < pSrc1->ncol; k++)
				{
					sum += pSrc1L[pSrc1->ncol * rowSrc1 + k] * pSrc2L[pSrc2->ncol * rowSrc2 + k];
				}
				pDstL[pDst->ncol * rowSrc1 + rowSrc2] = sum;
			}
		}
	}
	else
	{
		ResultL = MTX_SIZE_MISMATCH;
	}

	return ResultL;
}
enum mtxResultInfo mtx_chol_lower(Matrix_t *const pSrc)
{
	enum mtxResultInfo ResultL = MTX_OPERATION_OK;
	double *const pSrcL = pSrc->val;
	const uint8_t nrow = pSrc->nrow;
	const uint8_t ncol = pSrc->ncol;
	uint8_t col, row;
	int8_t tmp;
	double sum = 0;

	if (ncol == nrow)
	{
		const uint8_t mtxSize = nrow;

		for (col = 0; col < mtxSize; col++)
		{
			for (row = 0; row < mtxSize; row++)
			{
				sum = pSrcL[mtxSize * col + row];

				for (tmp = (int8_t)(col - 1); tmp >= 0; tmp--)
				{
					sum -= pSrcL[mtxSize * row + tmp] * pSrcL[mtxSize * col + tmp];
				}

				pSrcL[ncol * row + col] = (row == col) ? sqrt(sum) : (row > col) ? (sum / pSrcL[ncol * col + col]) : 0;

				if ((row == col) && (sum <= 0))
				{
					ResultL = MTX_NOT_POS_DEFINED;
				}
			}
		}
	}
	else
	{
		ResultL = MTX_NOT_SQUARE;
	}

	return ResultL;
}
enum mtxResultInfo mtx_chol_upper(Matrix_t *const pSrc)
{
	enum mtxResultInfo ResultL = MTX_OPERATION_OK;
	double *const pSrcL = pSrc->val;
	const uint8_t nrow = pSrc->nrow;
	const uint8_t ncol = pSrc->ncol;
	uint8_t col, row;
	int8_t tmp;
	double sum = 0;

	if (ncol == nrow)
	{
		for (row = 0; row < nrow; row++)
		{
			for (col = 0; col < ncol; col++)
			{
				sum = pSrcL[ncol * row + col];

				for (tmp = (int8_t)(row - 1); tmp >= 0; tmp--) // tmp could be calc negative
				{
					sum -= pSrcL[ncol * tmp + row] * pSrcL[ncol * tmp + col];
				}

				pSrcL[ncol * row + col] = (row == col) ? sqrt(sum) : (row < col) ? (sum / pSrcL[ncol * row + row]) : 0;

				if ((row == col) && (sum <= 0))
				{
					ResultL = MTX_NOT_POS_DEFINED;
				}
			}
		}
	}
	else
	{
		ResultL = MTX_NOT_SQUARE;
	}

	return ResultL;
}
/**
 * Upper cholesky decomposition variant 1
 */
enum mtxResultInfo mtx_chol1(double *A, double *L, uint8_t size)
{
	uint8_t Result = MTX_OPERATION_OK;
	uint8_t col, row;
	int8_t tmp;
	double sum = 0;

	for (row = 0; row < size; row++)
	{
		for (col = 0; col < size; col++)
		{
			sum = A[size * row + col];

			for (tmp = (int8_t)(row - 1); tmp >= 0; tmp--)
			{
				sum -= L[size * tmp + row] * L[size * tmp + col];
			}

			if (row == col)
			{
				if (sum > 0)
				{
					L[size * row + col] = sqrt(sum);
				}
				else
				{
					Result = MTX_NOT_POS_DEFINED;
				}
			}
			else if (row < col)
			{
				L[size * row + col] = sum / L[size * row + row]; // sum/Lii(diag)
			}
			else
			{
				L[size * row + col] = 0;
			}
		}
	}

	return Result;
}
/**
 * @brief Matrix inv  \f$A=A^{-1}\f$
 *
 * @param pDst At the begining should point to identity matrix!!
 * @param pSrc Square matrix
 */
enum mtxResultInfo mtx_inv(Matrix_t *const pSrc, Matrix_t *const pDst)
{
	enum mtxResultInfo Result = MTX_OPERATION_OK;
	const uint8_t nrow = pSrc->nrow;
	const uint8_t ncol = pSrc->ncol;
	uint8_t j, i;
	uint8_t k = 0;
	uint8_t l = 0;
	double s = 0;
	double t = 0;

	if (nrow == ncol)
	{
		for (j = 0; j < nrow; j++)
		{
			for (i = j; i < nrow; i++)
			{
				if (0 != pSrc->val[ncol * i + j])
				{
					for (k = 0; k < nrow; k++)
					{
						s = pSrc->val[ncol * j + k];
						pSrc->val[ncol * j + k] = pSrc->val[ncol * i + k];
						pSrc->val[ncol * i + k] = s;

						s = pDst->val[ncol * j + k];
						pDst->val[ncol * j + k] = pDst->val[ncol * i + k];
						pDst->val[ncol * i + k] = s;
					}

					t = 1 / pSrc->val[ncol * j + j];

					for (k = 0; k < nrow; k++)
					{
						pSrc->val[ncol * j + k] = t * pSrc->val[ncol * j + k];
						pDst->val[ncol * j + k] = t * pDst->val[ncol * j + k];
					}

					for (l = 0; l < nrow; l++)
					{
						if (l != j)
						{
							t = -pSrc->val[ncol * l + j];
							for (k = 0; k < nrow; k++)
							{
								pSrc->val[ncol * l + k] += t * pSrc->val[ncol * j + k];
								pDst->val[ncol * l + k] += t * pDst->val[ncol * j + k];
							}
						}
					}
				}
				break;
			}
			if (0 == pSrc->val[ncol * l + k])
			{
				Result = MTX_SINGULAR;
			}
		}
	}
	else
	{
		Result = MTX_SIZE_MISMATCH;
	}

	return Result;
}
enum mtxResultInfo mtx_add(Matrix_t *const pDst, Matrix_t const *const pSrc)
{
	uint8_t Result = MTX_OPERATION_OK;
	double *const pDstL = (double *)pDst->val;
	double const *const pSrcL = (double *)pSrc->val;
	uint16_t eIdx;

	if (pDst->ncol == pSrc->ncol && pDst->nrow == pSrc->nrow)
	{
		for (eIdx = 0; eIdx < pSrc->nelem; eIdx++)
		{
			pDstL[eIdx] += pSrcL[eIdx];
		}
	}
	else
	{
		Result = MTX_SIZE_MISMATCH;
	}

	return Result;
}
enum mtxResultInfo mtx_sub(Matrix_t *const pDst, Matrix_t const *const pSrc)
{
	uint8_t Result = MTX_OPERATION_OK;
	double *const pDstL = (double *)pDst->val;
	double const *const pSrcL = (double *)pSrc->val;
	uint16_t eIdx;

	if (pDst->ncol == pSrc->ncol && pDst->nrow == pSrc->nrow)
	{
		for (eIdx = 0; eIdx < pSrc->nelem; eIdx++)
		{
			pDstL[eIdx] -= pSrcL[eIdx];
		}
	}
	else
	{
		Result = MTX_SIZE_MISMATCH;
	}

	return Result;
}
enum mtxResultInfo mtx_mul_scalar(Matrix_t *const pSrc, const double scalar)
{
	enum mtxResultInfo Result = MTX_OPERATION_OK;
	double *const pDst = pSrc->val;
	uint16_t eIdx;

	for (eIdx = 0; eIdx < pSrc->nelem; eIdx++)
	{
		pDst[eIdx] *= scalar;
	}

	return Result;
}
enum mtxResultInfo mtx_sub_scalar(Matrix_t *const pSrc, const double scalar)
{
	enum mtxResultInfo Result = MTX_OPERATION_OK;
	double *const pDst = pSrc->val;
	uint16_t eIdx;

	for (eIdx = 0; eIdx < pSrc->nelem; eIdx++)
	{
		pDst[eIdx] -= scalar;
	}

	return Result;
}
enum mtxResultInfo mtx_add_scalar(Matrix_t *const pSrc, const double scalar)
{
	enum mtxResultInfo Result = MTX_OPERATION_OK;
	double *const pDst = pSrc->val;
	uint16_t eIdx;

	for (eIdx = 0; eIdx < pSrc->nelem; eIdx++)
	{
		pDst[eIdx] += scalar;
	}

	return Result;
}
enum mtxResultInfo mtx_cpy(Matrix_t *const pDst, Matrix_t const *const pSrc)
{
	enum mtxResultInfo Result = MTX_OPERATION_OK;
	double *const pDstL = pDst->val;
	double const *const pSrcL = pSrc->val;
	uint16_t eIdx;

	if (pDst->ncol == pSrc->ncol && pDst->nrow == pSrc->nrow)
	{
		for (eIdx = 0; eIdx < pSrc->nelem; eIdx++)
		{
			pDstL[eIdx] = pSrcL[eIdx];
		}
	}
	else
	{
		Result = MTX_SIZE_MISMATCH;
	}

	return Result;
}
enum mtxResultInfo mtx_identity(Matrix_t *const pSrc)
{
	enum mtxResultInfo Result = MTX_OPERATION_OK;
	double *const pDst = (double *)pSrc->val;
	const uint8_t nCol = pSrc->ncol;
	uint16_t eIdx;

	if (pSrc->nrow == nCol)
	{
		pDst[0] = 1;

		for (eIdx = 1; eIdx < pSrc->nelem; eIdx++)
		{
			const uint16_t cmpLeft = (uint16_t)(eIdx / nCol);

			pDst[eIdx] = eIdx < nCol ? 0 : cmpLeft == eIdx % (cmpLeft * nCol) ? 1 : 0;
		}
	}
	else
	{
		Result = MTX_SIZE_MISMATCH;
	}

	return Result;
}
enum mtxResultInfo mtx_zeros(Matrix_t *const pSrc)
{
	enum mtxResultInfo Result = MTX_OPERATION_OK;
	double *const pDst = (double *)pSrc->val;
	uint16_t eIdx;

	for (eIdx = 0; eIdx < pSrc->nelem; eIdx++)
	{
		pDst[eIdx] = 0;
	}

	return Result;
}
