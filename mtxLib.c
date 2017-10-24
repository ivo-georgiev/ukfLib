#include "mtxLib.h"

mtxResultInfo mtx_init_f64(tMatrix * A, float64 * pValue, uint8 nrow, uint8 ncol)
{
    A->val = pValue;
    A->ncol = ncol;
    A->nrow = nrow;

    return MTX_OPERATION_OK; 
}


mtxResultInfo mtx_diagsum_f64(tMatrix * pA, float64 * diagsum)
{
    mtxResultInfo Result = MTX_OPERATION_OK;
    const float64 * const pSrcL = (float64 *)pA->val;
    const uint8 nrow = pA->nrow;
    const uint8 ncol = pA->ncol;
    uint8 row,col;
    float64 sum=0;

    for(row=1;row<nrow;row++)//?bug this is for upper diag sum???
    {
        for(col=0;col<ncol;col++)
        {
            if(row > col)
            {
                sum += pSrcL[ncol*row+col];
            }
        }
    }
    *diagsum = sum;
    
    return Result;
}


//A=A'
mtxResultInfo mtx_transp_square_f64(tMatrix * pA)
{
    mtxResultInfo ResultL = MTX_OPERATION_OK;
    const uint8 nrow = pA->nrow;
    const uint8 ncol = pA->ncol;
    float64 * pSrcL = (float64 *)pA->val;
    uint8 row,col;
    float64 temp;
    
    if(nrow == ncol)
    {
        for(row=0;row<nrow;row++)
        {
            for(col=0;col<ncol;col++) 
            {
                if(row != col && row<col)
                {
                    temp = pSrcL[nrow*row+col];
                    pSrcL[ncol*row+col] = pSrcL[ncol*col+row];
                    pSrcL[ncol*col+row] = temp;     
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

//
mtxResultInfo mtx_transp_dest_f64(tMatrix * pA,tMatrix * pB)
{
    mtxResultInfo ResultL = MTX_OPERATION_OK;
    float64 * pSrcL = (float64 *)pA->val;
    float64 * pDstL = (float64 *)pB->val;
    const uint8 nRowSrcL = pA->nrow;
    const uint8 nColSrcL = pA->ncol;
    const uint8 nRowDstL = pB->nrow;
    const uint8 nColDstL = pB->ncol;
    uint8 row,col;
    
    if(nRowSrcL == nColDstL || nColSrcL == nRowDstL)
    {
        for(row=0;row<nRowDstL;row++)
        {
            for(col=0;col<nColDstL;col++)
            {
                pDstL[nColDstL*row + col] = pSrcL[nColSrcL*col + row];          
            }
        }
    }
    else
    {
        ResultL = MTX_SIZE_MISMATCH;   
    }

return ResultL;
}

//C=A*B
mtxResultInfo mtx_mul_f64(tMatrix * pA, tMatrix * pB, tMatrix * pC)
{
    mtxResultInfo ResultL = MTX_OPERATION_OK;
    float64 * const pSrc1L = (float64 *)pA->val;
    float64 * const pSrc2L = (float64 *)pB->val;
    float64 * const pDstL = (float64 *)pC->val;
    //const uint8 nrow = pA->nrow;
    //const uint8 ncol = pA->ncol;
    uint8 row,col,k;
    float64 sum;

    if(pA->ncol == pB->nrow)//?
    {        
        for(row=0;row<pA->nrow;row++)
        {
            for(col=0;col<pB->ncol;col++)  
            {
                sum = 0;
                for(k=0;k<pA->ncol;k++)
                {
                    sum += pSrc1L[pA->ncol*row+k] * pSrc2L[pB->ncol*k+col];
                }
                pDstL[pC->ncol*row+col] = sum;
            }
        }
    }
    else
    {
        ResultL = MTX_SIZE_MISMATCH;
    }
    
    return ResultL;
}
//A=chol_lower(A)
mtxResultInfo mtx_chol_lower_f64(tMatrix * pA)
{
    mtxResultInfo ResultL = MTX_OPERATION_OK;
    float64 * const pSrcL = pA->val;
    const uint8 nrow = pA->nrow;
    const uint8 ncol = pA->ncol;
    uint8 col,row;
    sint8 tmp;// 
    float64 sum=0;
    
    if(ncol == nrow)
    {   
        const uint8 mtxSize = nrow;

        for(col=0;col<mtxSize;col++)
        {
            for(row=0;row<mtxSize;row++)
            {
                sum = pSrcL[mtxSize*col+row];
                
                for(tmp = col-1;tmp>=0;tmp--)
                {
                    sum -= pSrcL[mtxSize*row+tmp] * pSrcL[mtxSize*col+tmp];
                }
                
                pSrcL[ncol*row + col] = (row==col) ? sqrt(sum) : (row > col) ? (sum / pSrcL[ncol*col+col]) : 0;
                
                
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
//A=chol_upper(A)
mtxResultInfo mtx_chol_upper_f64(tMatrix * pA)
{
    mtxResultInfo ResultL = MTX_OPERATION_OK;
    float64 * const pSrcL = pA->val;
    const uint8 nrow = pA->nrow;
    const uint8 ncol = pA->ncol;
    uint8 col,row;
    sint8 tmp;// 
    float64 sum=0;
    
    if(ncol == nrow)
    {       
        for(row=0;row<nrow;row++)
        {
            for(col=0;col<ncol;col++)
            {
                sum = pSrcL[ncol*row + col];
                
                for(tmp = row-1;tmp>=0;tmp--)// tmp could be calc negative
                {
                    sum -= pSrcL[ncol*tmp+row] * pSrcL[ncol*tmp+col];
                }
                
                pSrcL[ncol*row + col] = (row==col) ? sqrt(sum) : (row < col) ? (sum / pSrcL[ncol*row+row]) : 0;
                
                
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
mtxResultInfo mtx_chol1_f64(float64* A, float64* L,uint8 size)
{
    uint8 Result = MTX_OPERATION_OK;
    uint8 col,row;
    sint8 tmp;
    float64 sum=0;
    
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

int mtx_sum_updiag_f64(int*A ,int m, int n)
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

//* pA - source  pI- destination
mtxResultInfo mtx_inv_f64(tMatrix * pA, tMatrix * pI)
{
    mtxResultInfo ResultL = MTX_OPERATION_OK;
    const uint8 nrow = pA->nrow;
    const uint8 ncol = pA->ncol;
    uint8 j,i,k,l;
    float64 s=0;
    float64 t=0;
    
    for(j = 0;j<nrow;j++)
    {
        for(i = j; i<nrow; i++)
        {
            if(0 != pA->val[ncol*i+j])
            {
                for(k = 0;k<nrow;k++)
                {
                    s = pA->val[ncol*j+k];
                    pA->val[ncol*j+k] = pA->val[ncol*i+k];
                    pA->val[ncol*i+k] = s;
                    
                    s = pI->val[ncol*j+k];
                    pI->val[ncol*j+k] = pI->val[ncol*i+k];
                    pI->val[ncol*i+k] = s;
                }
                
                t = 1 / pA->val[ncol*j+j];
                
                for(k=0;k<nrow;k++)
                {
                    pA->val[ncol*j+k] = t * pA->val[ncol*j+k];
                    pI->val[ncol*j+k] = t * pI->val[ncol*j+k];
                }

                for(l=0;l<nrow;l++)
                {
                    if(l != j)
                    {
                        t = -pA->val[ncol*l+j];
                        for(k=0;k<nrow;k++)
                        {
                            pA->val[ncol*l+k] += t *  pA->val[ncol*j+k];
                            pI->val[ncol*l+k] += t *  pI->val[ncol*j+k];
                        }
                    }
                }
            }
            break;
        }
        if(0 == pA->val[ncol*l+k])
        {
            ResultL = MTX_SINGULAR;
        }
    }
    
    return ResultL;
    
}

//A=A+B
mtxResultInfo mtx_add_f64(tMatrix * pA,tMatrix * pB)
{
    mtxResultInfo Result = MTX_OPERATION_OK;
    float64 * pDest= (float64 *)pA->val;
    float64 * pSrc1= (float64 *)pB->val;

    const uint8 nRow = pA->nrow;
    const uint8 nCol = pA->ncol;
    uint8 col,row;

    // float64 sum=0;
    
    for(row=0;row<nRow;row++)
    {
        for(col=0;col<nCol;col++)
        {
            pDest[nCol*row+col] += pSrc1[nCol*row+col];
        }
    }

    return Result;
}
//A=A-B
mtxResultInfo mtx_subtract_f64(tMatrix * pA,tMatrix * pB)
{
    uint8 Result = MTX_OPERATION_OK;
    float64 * pDest= (float64 *)pA->val;
    float64 * pSrc1= (float64 *)pB->val;    
    const uint8 nRow = pA->nrow;
    const uint8 nCol = pA->ncol;
    uint8 col,row;
  
    
    for(row=0;row<nRow;row++)
    {
        for(col=0;col<nCol;col++)
        {
            pDest[nCol*row+col] -= pSrc1[nCol*row+col];
        }
    }
    
    return Result;
}

//A=k*A=A*k
mtxResultInfo mtx_mul_scalar_f64(tMatrix * pA,float64 scalar)
{
    mtxResultInfo Result = MTX_OPERATION_OK;
    float64 * pDest= pA->val;    
    const uint8 nRow = pA->nrow;
    const uint8 nCol = pA->ncol;
    uint8 col,row;
    
    
    for(row=0;row<nRow;row++)
    {
        for(col=0;col<nCol;col++)
        {
            pDest[nCol*row+col] *= scalar;
        }
    }
    
    return Result;
}

mtxResultInfo mtx_subtract_scalar_f64(tMatrix * pA,float64 scalar)
{
    mtxResultInfo Result = MTX_OPERATION_OK;
    float64 * pDest= pA->val;    
    const uint8 nRow = pA->nrow;
    const uint8 nCol = pA->ncol;
    uint8 col,row;
    
    
    for(row=0;row<nRow;row++)
    {
        for(col=0;col<nCol;col++)
        {
            pDest[nCol*row+col] -= scalar;
        }
    }
    
    return Result;
}

mtxResultInfo mtx_add_scalar_f64(tMatrix * pA,float64 scalar)
{
    mtxResultInfo Result = MTX_OPERATION_OK;
    float64 * pDest= pA->val;    
    const uint8 nRow = pA->nrow;
    const uint8 nCol = pA->ncol;
    uint8 col,row;
    
    
    for(row=0;row<nRow;row++)
    {
        for(col=0;col<nCol;col++)
        {
            pDest[nCol*row+col] += scalar;
        }
    }
    
    return Result;
}

mtxResultInfo mtx_cpy_f64(tMatrix * const pDestP,tMatrix const * const pSrcP)
{
    mtxResultInfo Result = MTX_OPERATION_OK;
    float64 * const pDestL = pDestP->val;
    float64 const * const pSrcL= pSrcP->val;
    const uint8 nRow = pDestP->nrow;
    const uint8 nCol = pSrcP->ncol;
    uint8 row,col;

    if(pDestP->ncol == pSrcP->ncol && pDestP->nrow == pSrcP->nrow)
    {
        for(row=0;row<nRow;row++)
        {
            for(col=0;col<nCol;col++)
            {
                pDestL[nCol*row+col] = pSrcL[nCol*row+col];
            }
        }
        
        //memcpy(pDest,pSrc,sizeof(pSrc));
    }
    else
    {
        Result = MTX_SIZE_MISMATCH;

    }

    return Result;

}

mtxResultInfo mtx_identity_f64(tMatrix * pI)
{
    mtxResultInfo Result = MTX_OPERATION_OK;
    float64 * pDest= (float64 *)pI->val;    
    const uint8 nRow = pI->nrow;
    const uint8 nCol = pI->ncol;
    uint8 col,row;
    
    if(nRow == nCol)
    {        
        for(row=0;row<nRow;row++)
        {
            for(col=0;col<nCol;col++)
            {
                pDest[nCol*row+col] = (row==col) ?  1 : 0;             
            }
        }
    }
    else
    {
        Result = MTX_SIZE_MISMATCH;
    }

    return Result;
}

mtxResultInfo mtx_zeros_f64(tMatrix * pZ)
{
    mtxResultInfo Result = MTX_OPERATION_OK;
    float64 * pDest= (float64 *)pZ->val;    
    const uint8 nRow = pZ->nrow;
    const uint8 nCol = pZ->ncol;
    uint8 col,row;    
    
    for(row=0;row<nRow;row++)
    {
        for(col=0;col<nCol;col++)
        {
            pDest[nCol*row+col] = 0;             
        }
    }
    
    return Result;
}
