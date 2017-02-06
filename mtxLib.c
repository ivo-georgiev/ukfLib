#include<math.h>
#include "mtxLib.h"


mtxResultInfo mtxLib_Transp_dp(sMatrixType A)
{
    mtxResultInfo ResultL = MTX_OPERATION_OK;
    int i,j;
    double temp;

    for(i=0;i<A.row;i++)
    {
        for(j=0;j<A.col;j++)
        {
            if(i != j && i<j)
            {
                temp = *(A.val+i*A.col+j);
                *(A.val+i*A.col+j) = *(A.val+j*A.col+i);
                *(A.val+j*A.col+i) = temp;
            }
        }
    }
    return ResultL;
}


mtxResultInfo mtxLib_Mul_dp(sMatrixType A, sMatrixType B, sMatrixType C)
{
    mtxResultInfo ResultL = MTX_OPERATION_OK;
    int i,j,k;
    double sum;
    
    if(A.col != B.row)
    {
        ResultL = MTX_OPERATION_ERROR;
    }
    for(i=0;i<A.row;i++)//A.row= 2
    {
        for(j=0;j<B.col;j++)//B.col  
        {
            sum = 0;
            for(k=0;k<B.row;k++)//B.col = 2
            {
                //sum = sum+ A[i][k] * B[k][j];
                sum += *(A.val+i*A.col+k) * *(B.val+k*B.col+j);
            }
            *(C.val+i*C.col+j) = sum;//C[i][j] 
        }
    }
    return ResultL;
}


void mtxLib_Cholesky_LL_dp(double* mtxA,int sizeA)
{
    int colJ, rowI,tmpK;
    double sumL=0;
    
    for(rowI=0;rowI<sizeA;rowI++)
    {
        for(colJ=0;colJ<sizeA;colJ++)
        {
            sumL = *(mtxA+sizeA*rowI+colJ);
            
            for(tmpK = rowI-1;tmpK>=0;tmpK--)
            {
                sumL -= *(mtxA+sizeA*tmpK+rowI) *  *(mtxA+sizeA*tmpK+colJ);
            }
            
            *(mtxA+sizeA*rowI+colJ) = (rowI==colJ) ? sqrt(sumL): (sumL/ *(mtxA+sizeA*rowI + rowI));
        }
    }
}

//cholesky Decomposition Upper variant 1
int mtxLib_Cholesky1_LL_dp(double* mtxA, double* mtxL,int mtxSize)
{
    int colIdxL, rowIdxL,tempIdxL;
    int ReturnL = 0;
    double sumL=0;
    
    for(rowIdxL=0;rowIdxL<mtxSize;rowIdxL++)
    {
        for(colIdxL=0;colIdxL<mtxSize;colIdxL++)
        {
            sumL = *(mtxA+mtxSize*rowIdxL+colIdxL);
            
            for(tempIdxL = rowIdxL-1;tempIdxL>=0;tempIdxL--)
            {
                sumL -= *(mtxL+mtxSize*tempIdxL+rowIdxL) *  *(mtxL+mtxSize*tempIdxL+colIdxL);
            }
            
            if(rowIdxL==colIdxL)
            {
                if(sumL>0)
                {
                    *(mtxL+mtxSize*rowIdxL + colIdxL) = sqrt(sumL);
                }
                else
                {
                    //matrix is not positive definite
                    ReturnL = 1;
                }
            }
            else if(rowIdxL < colIdxL)
            {
                *(mtxL+mtxSize*rowIdxL + colIdxL) = sumL/ *(mtxL+mtxSize*rowIdxL + rowIdxL);
            }
            else
            {
                *(mtxL+mtxSize*rowIdxL + colIdxL) = 0;
            }
        }
    }
    
    return ReturnL;
}