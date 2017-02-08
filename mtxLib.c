#include<math.h>
#include "mtxLib.h"

mtxResultInfo mtx_diagsum_dp(sMatrixType A, double * DiagSum)
{
    mtxResultInfo ResultL = MTX_OPERATION_OK;
    int i,j;
    double sum=0;

    for(i=1;i<A.nrow;i++)
    {
        for(j=0;j<A.ncol;j++)
        {
            if(i > j)
            {
                sum += *(A.val+i*A.ncol+j);
            }
        }
    }
    *DiagSum = sum;
    
    return ResultL;
}



//A=A'
mtxResultInfo mtxLib_Transp_dp(sMatrixType A)
{
    mtxResultInfo ResultL = MTX_OPERATION_OK;
    int i,j;
    double temp;

    for(i=0;i<A.nrow;i++)
    {
        for(j=0;j<A.ncol;j++)
        {
            if(i != j && i<j)
            {
                temp = *(A.val+i*A.ncol+j);
                *(A.val+i*A.ncol+j) = *(A.val+j*A.ncol+i);
                *(A.val+j*A.ncol+i) = temp;
            }
        }
    }
    return ResultL;
}

//C=A*B
mtxResultInfo mtxLib_Mul_dp(sMatrixType A, sMatrixType B, sMatrixType C)
{
    mtxResultInfo ResultL = MTX_OPERATION_OK;
    int row,col,k;
    double sum;
    
    if(A.ncol != B.nrow)
    {
        ResultL = MTX_OPERATION_ERROR;
    }
    for(row=0;row<A.nrow;row++)//A.nrow= 2
    {
        for(col=0;col<B.ncol;col++)//B.col  
        {
            sum = 0;
            for(k=0;k<B.nrow;k++)//B.col = 2
            {
                //sum = sum+ A[i][k] * B[k][col];
                sum += *(A.val+row*A.ncol+k) * *(B.val+k*B.ncol+col);//*((int *)y + 2 * NUMBER_OF_COLUMNS + 2); // Right!
            }
            *(C.val+row*C.ncol+col) = sum;//C[i][col] 
        }
    }
    return ResultL;
}

//A=chol(A)
//http://rosettacode.org/wiki/Cholesky_decomposition
mtxResultInfo mtxLib_Cholesky_LL_dp(sMatrixType A)
{
    mtxResultInfo ResultL = MTX_OPERATION_OK;
    int col, row,tmp;
    double sumL=0;
    
    if(A.ncol == A.nrow)
    {       
        for(row=0;row<A.nrow;row++)
        {
            for(col=0;col<A.ncol;col++)//col<A.ncol
            {
                sumL = A.val[A.ncol*row + col];
                
                for(tmp = row-1;tmp>=0;tmp--)
                {
                    sumL -= A.val[A.ncol*tmp+row] *  A.val[A.ncol*tmp+col];
                }
                
                A.val[A.ncol*row + col] = (row==col) ? sqrt(sumL): (sumL/A.val[A.ncol*row+row]);
            }
        }
    }
    else
    {
        ResultL = MTX_OPERATION_ERROR;//not square
    }

    return ResultL;
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