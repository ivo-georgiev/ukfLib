#include<math.h>
#include "mtxLib.h"


mtxResultInfo mtx_init_dp(sMatrixType* A, double * M, int nrow, int ncol)
{
    A->val = M;
    A->ncol = ncol;
    A->nrow = nrow;

    return MTX_OPERATION_OK; 
}


mtxResultInfo mtx_diagsum_dp(sMatrixType A, double * diagsum)
{
    mtxResultInfo Result = MTX_OPERATION_OK;
    int row,col;
    double sum=0;

    for(row=1;row<A.nrow;row++)
    {
        for(col=0;col<A.ncol;col++)
        {
            if(row > col)
            {
                sum += A.val[A.ncol*row+col];
            }
        }
    }
    *diagsum = sum;
    
    return Result;
}


//A=A'
mtxResultInfo mtx_transp_dp(sMatrixType A)
{
    mtxResultInfo ResultL = MTX_OPERATION_OK;
    int row,col;
    double temp;

    for(row=0;row<A.nrow;row++)
    {
        for(col=0;col<A.ncol;col++)
        {
            if(row != col && row<col)
            {
                temp = A.val[A.ncol*row+col];
                A.val[A.ncol*row+col] = A.val[A.ncol*col+row];
                A.val[A.ncol*col+row] = temp;
            }
        }
    }
    return ResultL;
}

//C=A*B
mtxResultInfo mtx_mul_dp(sMatrixType A, sMatrixType B, sMatrixType C)
{
    mtxResultInfo ResultL = MTX_OPERATION_OK;
    int row,col,k;
    double sum;
    
    if(A.ncol != B.nrow)
    {
        ResultL = MTX_DIMENSION_ERR;
    }
    for(row=0;row<A.nrow;row++)
    {
        for(col=0;col<B.ncol;col++)  
        {
            sum = 0;
            for(k=0;k<B.nrow;k++)
            {
                sum += A.val[A.ncol*row+k] * B.val[B.ncol*k+col];
            }
            C.val[C.ncol*row+col] = sum;
        }
    }
    return ResultL;
}

//A=chol(A)
//http://rosettacode.org/wiki/Cholesky_decomposition
mtxResultInfo mtx_chol_dp(sMatrixType A)
{
    mtxResultInfo ResultL = MTX_OPERATION_OK;
    int col,row,tmp;
    double sum=0;
    
    if(A.ncol == A.nrow)
    {       
        for(row=0;row<A.nrow;row++)
        {
            for(col=0;col<A.ncol;col++)
            {
                sum = A.val[A.ncol*row + col];
                
                for(tmp = row-1;tmp>=0;tmp--)
                {
                    sum -= A.val[A.ncol*tmp+row] * A.val[A.ncol*tmp+col];
                }
                
                A.val[A.ncol*row + col] = (row==col) ? sqrt(sum) : (sum / A.val[A.ncol*row+row]);
                
                
                if((row==col) && (sum<=0))
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

//cholesky Decomposition Upper variant 1
int mtx_chol1_dp(double* A, double* L,int size)
{
    int Result = MTX_OPERATION_OK;
    int col,row,tmp;
    double sum=0;
    
    for(row=0;row<size;row++)
    {
        for(col=0;col<size;col++)
        {
            sum = A[size*row + col];
            
            for(tmp = row-1;tmp>=0;tmp--)
            {
                sum -= L[size*tmp + row] * L[size*tmp + col];
            }
            
            if(row==col)
            {
                if(sum>0)
                {
                    L[size*row + col] = sqrt(sum);
                }
                else
                {
                    Result = MTX_NOT_POS_DEFINED;
                }
            }
            else if(row < col)
            {
                L[size*row + col] = sum/ L[size*row + row];//sum/Lii(diag)
            }
            else
            {
                L[size*row + col] = 0;
            }
        }
    }
    
    return Result;
}

int mtx_sum_updiag_dp(int*A ,int m, int n)
{
    int i,j,sum=0;
    for(i=0;i<m;i++)
    {
        for(j=1;j<n;j++)
        {
            if(i < j)
            {
                sum += *(A+i*n+j);
            }
        }
    }
    return sum;
}