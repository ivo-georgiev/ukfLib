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

/******************************************************************************************************************************************************************************************************\
 ***  FUNCTION:
 ***     
 *** 
 ***  DESCRIPTION:
 ***
 ***            
 ***  PARAMETERS:
 ***      Type               Name              Range              Description
 ***      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 ***      
 ***  RETURNS:
 ***      
 ***  SETTINGS:
 ***
\******************************************************************************************************************************************************************************************************/
mtxResultInfo mtx_init_bool(tMatrixBool * const pSrc, boolean * const pValue, const uint8 nrow, const uint8 ncol,const uint16 nelem)
{
    pSrc->val = pValue;
    pSrc->ncol = ncol;
    pSrc->nrow = nrow;
    pSrc->nelem = nelem;
    return MTX_OPERATION_OK; 
}
/******************************************************************************************************************************************************************************************************\
 ***  FUNCTION:
 ***     
 *** 
 ***  DESCRIPTION:
 ***
 ***            
 ***  PARAMETERS:
 ***      Type               Name              Range              Description
 ***      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 ***      
 ***  RETURNS:
 ***      
 ***  SETTINGS:
 ***
\******************************************************************************************************************************************************************************************************/
mtxResultInfo mtx_init_f64(tMatrix * const pSrc, float64 * const pValue, const uint8 nrow, const uint8 ncol,const uint16 nelem)
{
    pSrc->val = pValue;
    pSrc->ncol = ncol;
    pSrc->nrow = nrow;
    pSrc->nelem = nelem;
    return MTX_OPERATION_OK; 
}
/******************************************************************************************************************************************************************************************************\
 ***  FUNCTION:
 ***     
 *** 
 ***  DESCRIPTION: For square matrix only
 ***
 ***            
 ***  PARAMETERS:
 ***      Type               Name              Range              Description
 ***      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 ***      
 ***  RETURNS:
 ***      
 ***  SETTINGS:
 ***
\******************************************************************************************************************************************************************************************************/
mtxResultInfo mtx_diagsum_f64(tMatrix * pSrc, float64 * diagsum)
{
    mtxResultInfo Result = MTX_OPERATION_OK;
    float64 const * const pSrcL = (float64 *)pSrc->val;
    const uint8 ncol = pSrc->ncol;
    uint16 eIdx;
    float64 sum = pSrcL[0];
    
    if(pSrc->nrow == ncol)
    {
        for(eIdx=1;eIdx < pSrc->nelem;eIdx++)
        {
            const uint16 cmpLeft = (uint16)(eIdx / ncol);
            
            sum += eIdx < ncol ? 0 : cmpLeft == eIdx % (cmpLeft*ncol) ?  pSrcL[eIdx] : 0;   
        }
    }
    else
    {
        Result = MTX_SIZE_MISMATCH;
    }
 

    *diagsum = sum;
    
    return Result;
}
/******************************************************************************************************************************************************************************************************\
 ***  FUNCTION:
 ***     
 *** 
 ***  DESCRIPTION: A=A'
 ***
 ***            
 ***  PARAMETERS:
 ***      Type               Name              Range              Description
 ***      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 ***      
 ***  RETURNS:
 ***      
 ***  SETTINGS:
 ***
\******************************************************************************************************************************************************************************************************/
mtxResultInfo mtx_transp_square_f64(tMatrix * const pSrc)
{
    mtxResultInfo ResultL = MTX_OPERATION_OK;
    const uint8 nrow = pSrc->nrow;
    const uint8 ncol = pSrc->ncol;
    float64 * const pSrcL = (float64 *)pSrc->val;
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
/******************************************************************************************************************************************************************************************************\
 ***  FUNCTION:
 ***     
 *** 
 ***  DESCRIPTION:
 ***
 ***            
 ***  PARAMETERS:
 ***      Type               Name              Range              Description
 ***      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 ***      
 ***  RETURNS:
 ***      
 ***  SETTINGS:
 ***
\******************************************************************************************************************************************************************************************************/
mtxResultInfo mtx_transp_dest_f64(tMatrix const * const pSrc,tMatrix * const pDst)
{
    mtxResultInfo ResultL = MTX_OPERATION_OK;
    float64 const * const pSrcL = (float64 *)pSrc->val;
    float64 * const pDstL = (float64 *)pDst->val;
    const uint8 nRowSrcL = pSrc->nrow;
    const uint8 nColSrcL = pSrc->ncol;
    const uint8 nRowDstL = pDst->nrow;
    const uint8 nColDstL = pDst->ncol;
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
/******************************************************************************************************************************************************************************************************\
 ***  FUNCTION:
 ***     
 *** 
 ***  DESCRIPTION: C=A*B
 ***
 ***            
 ***  PARAMETERS:
 ***      Type               Name              Range              Description
 ***      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 ***      
 ***  RETURNS:
 ***      
 ***  SETTINGS:
 ***
\******************************************************************************************************************************************************************************************************/
mtxResultInfo mtx_mul_f64(tMatrix const * const pSrc1, tMatrix const * const pSrc2, tMatrix * const pDst)
{
    mtxResultInfo ResultL = MTX_OPERATION_OK;
    float64 const * const pSrc1L = (float64 *)pSrc1->val;
    float64 const * const pSrc2L = (float64 *)pSrc2->val;
    float64 * const pDstL = (float64 *)pDst->val;
    uint8 row,col,k;
    float64 sum;

    if(pSrc1->ncol == pSrc2->nrow)
    {        
        for(row=0;row<pSrc1->nrow;row++)
        {
            for(col=0;col<pSrc2->ncol;col++)  
            {
                sum = 0;
                for(k=0;k<pSrc1->ncol;k++)
                {
                    sum += pSrc1L[pSrc1->ncol*row+k] * pSrc2L[pSrc2->ncol*k+col];
                }
                pDstL[pDst->ncol*row+col] = sum;
            }
        }
    }
    else
    {
        ResultL = MTX_SIZE_MISMATCH;
    }
    
    return ResultL;
}
/******************************************************************************************************************************************************************************************************\
 ***  FUNCTION: 
 ***     
 *** 
 ***  DESCRIPTION: Special multiplication Dst=Src1*Src2'
 ***  Special function for multiplication of Src1 matrix with transpose image of Src2 matrix. Be sure that array size are suitable for multiplication after Src2 transpose  
 ***            
 ***  PARAMETERS:
 ***      Type               Name              Range              Description
 ***      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 ***      
 ***  RETURNS:
 ***      
 ***  SETTINGS:
 ***
\******************************************************************************************************************************************************************************************************/
mtxResultInfo mtx_mul_src2tr_f64(tMatrix const * const pSrc1, tMatrix const * const pSrc2, tMatrix * const pDst)
{
    mtxResultInfo ResultL = MTX_OPERATION_OK;
    float64 const * const pSrc1L = (float64 *)pSrc1->val;
    float64 const * const pSrc2L = (float64 *)pSrc2->val;
    float64 * const pDstL = (float64 *)pDst->val;
    uint8 rowSrc1,rowSrc2,k;
    float64 sum;

    if(pSrc1->ncol == pSrc2->ncol)
    {        
        for(rowSrc1=0;rowSrc1<pSrc1->nrow;rowSrc1++)
        {
            for(rowSrc2=0;rowSrc2<pSrc2->nrow;rowSrc2++)  
            {
                sum = 0;
                for(k=0;k<pSrc1->ncol;k++)
                {
                    sum += pSrc1L[pSrc1->ncol*rowSrc1 + k ] * pSrc2L[pSrc2->ncol*rowSrc2 + k];
                }
                pDstL[pDst->ncol*rowSrc1 + rowSrc2] = sum;
            }
        }
    }
    else
    {
        ResultL = MTX_SIZE_MISMATCH;
    }
    
    return ResultL;
}
/******************************************************************************************************************************************************************************************************\
 ***  FUNCTION:
 ***     
 *** 
 ***  DESCRIPTION:
 ***
 ***            
 ***  PARAMETERS:
 ***      Type               Name              Range              Description
 ***      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 ***      
 ***  RETURNS:
 ***      
 ***  SETTINGS:
 ***
\******************************************************************************************************************************************************************************************************/
mtxResultInfo mtx_chol_lower_f64(tMatrix * const pSrc)
{
    mtxResultInfo ResultL = MTX_OPERATION_OK;
    float64 * const pSrcL = pSrc->val;
    const uint8 nrow = pSrc->nrow;
    const uint8 ncol = pSrc->ncol;
    uint8 col,row;
    sint8 tmp; 
    float64 sum=0;
    
    if(ncol == nrow)
    {   
        const uint8 mtxSize = nrow;

        for(col=0;col<mtxSize;col++)
        {
            for(row=0;row<mtxSize;row++)
            {
                sum = pSrcL[mtxSize*col+row];
                
                for(tmp = (sint8)(col-1);tmp>=0;tmp--)
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
/******************************************************************************************************************************************************************************************************\
 ***  FUNCTION:
 ***     
 *** 
 ***  DESCRIPTION:
 ***
 ***            
 ***  PARAMETERS:
 ***      Type               Name              Range              Description
 ***      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 ***      
 ***  RETURNS:
 ***      
 ***  SETTINGS:
 ***
\******************************************************************************************************************************************************************************************************/
mtxResultInfo mtx_chol_upper_f64(tMatrix * const pSrc)
{
    mtxResultInfo ResultL = MTX_OPERATION_OK;
    float64 * const pSrcL = pSrc->val;
    const uint8 nrow = pSrc->nrow;
    const uint8 ncol = pSrc->ncol;
    uint8 col,row;
    sint8 tmp;
    float64 sum=0;
    
    if(ncol == nrow)
    {       
        for(row=0;row<nrow;row++)
        {
            for(col=0;col<ncol;col++)
            {
                sum = pSrcL[ncol*row + col];
                
                for(tmp = (sint8)(row-1);tmp>=0;tmp--)// tmp could be calc negative
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
/******************************************************************************************************************************************************************************************************\
 ***  FUNCTION:
 ***     
 *** 
 ***  DESCRIPTION: Upper cholesky decomposition variant 1
 ***
 ***            
 ***  PARAMETERS:
 ***      Type               Name              Range              Description
 ***      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 ***      
 ***  RETURNS:
 ***      
 ***  SETTINGS:
 ***
\******************************************************************************************************************************************************************************************************/
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
            
            for(tmp = (sint8)(row-1);tmp>=0;tmp--)
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
/******************************************************************************************************************************************************************************************************\
***  FUNCTION:
***     
 *** 
 ***  DESCRIPTION:
 ***  
 ***            
 ***  PARAMETERS:
 ***      Type               Name              Range              Description
 ***      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 ***      tMatrix * const    pDst                                - At the begining should point to identity matrix!!
 ***      tMatrix * const    pSrc                                - Square matrix
 ***  RETURNS:
 ***      
 ***  SETTINGS:
 ***
\******************************************************************************************************************************************************************************************************/
mtxResultInfo mtx_inv_f64(tMatrix * const pSrc, tMatrix * const pDst)
{
    mtxResultInfo Result = MTX_OPERATION_OK;
    const uint8 nrow = pSrc->nrow;
    const uint8 ncol = pSrc->ncol;
    uint8 j,i;
    uint8 k = 0;
    uint8 l = 0;
    float64 s=0;
    float64 t=0;
    
    if(nrow == ncol)
    {     
        for(j = 0;j<nrow;j++)
        {
            for(i = j; i<nrow; i++)
            {
                if(0 != pSrc->val[ncol*i+j])
                {
                    for(k = 0;k<nrow;k++)
                    {
                        s = pSrc->val[ncol*j+k];
                        pSrc->val[ncol*j+k] = pSrc->val[ncol*i+k];
                        pSrc->val[ncol*i+k] = s;
                        
                        s = pDst->val[ncol*j+k];
                        pDst->val[ncol*j+k] = pDst->val[ncol*i+k];
                        pDst->val[ncol*i+k] = s;
                    }
                    
                    t = 1 / pSrc->val[ncol*j+j];
                    
                    for(k=0;k<nrow;k++)
                    {
                        pSrc->val[ncol*j+k] = t * pSrc->val[ncol*j+k];
                        pDst->val[ncol*j+k] = t * pDst->val[ncol*j+k];
                    }
                    
                    for(l=0;l<nrow;l++)
                    {
                        if(l != j)
                        {
                            t = -pSrc->val[ncol*l+j];
                            for(k=0;k<nrow;k++)
                            {
                                pSrc->val[ncol*l+k] += t *  pSrc->val[ncol*j+k];
                                pDst->val[ncol*l+k] += t *  pDst->val[ncol*j+k];
                            }
                        }
                    }
                }
                break;
            }
            if(0 == pSrc->val[ncol*l+k])
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
/******************************************************************************************************************************************************************************************************\
 ***  FUNCTION:
 ***     
 *** 
 ***  DESCRIPTION:
 ***
 ***            
 ***  PARAMETERS:
 ***      Type               Name              Range              Description
 ***      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 ***      
 ***  RETURNS:
 ***      
 ***  SETTINGS:
 ***
\******************************************************************************************************************************************************************************************************/
mtxResultInfo mtx_add_f64(tMatrix * const pDst,tMatrix const * const pSrc)
{
    uint8 Result = MTX_OPERATION_OK;
    float64 * const pDstL= (float64 *)pDst->val;
    float64 const * const pSrcL= (float64 *)pSrc->val;    
    uint16 eIdx;
  
    if(pDst->ncol == pSrc->ncol && pDst->nrow == pSrc->nrow)
    {
        for(eIdx=0;eIdx<pSrc->nelem;eIdx++)
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
/******************************************************************************************************************************************************************************************************\
 ***  FUNCTION:
 ***     
 *** 
 ***  DESCRIPTION:
 ***
 ***            
 ***  PARAMETERS:
 ***      Type               Name              Range              Description
 ***      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 ***      
 ***  RETURNS:
 ***      
 ***  SETTINGS:
 ***
\******************************************************************************************************************************************************************************************************/
mtxResultInfo mtx_sub_f64(tMatrix * const pDst,tMatrix const * const pSrc)
{
    uint8 Result = MTX_OPERATION_OK;
    float64 * const pDstL= (float64 *)pDst->val;
    float64 const * const pSrcL= (float64 *)pSrc->val;    
    uint16 eIdx;
  
    if(pDst->ncol == pSrc->ncol && pDst->nrow == pSrc->nrow)
    {
        for(eIdx=0;eIdx<pSrc->nelem;eIdx++)
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
/******************************************************************************************************************************************************************************************************\
 ***  FUNCTION:
 ***     
 *** 
 ***  DESCRIPTION:
 ***
 ***            
 ***  PARAMETERS:
 ***      Type               Name              Range              Description
 ***      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 ***      
 ***  RETURNS:
 ***      
 ***  SETTINGS:
 ***
\******************************************************************************************************************************************************************************************************/
mtxResultInfo mtx_mul_scalar_f64(tMatrix * const pSrc,const float64 scalar)
{
    mtxResultInfo Result = MTX_OPERATION_OK;
    float64 * const pDst= pSrc->val;    
    uint16 eIdx;
    
    for(eIdx=0;eIdx<pSrc->nelem;eIdx++)
    {
        pDst[eIdx] *= scalar;
    }
    
    return Result;
}
/******************************************************************************************************************************************************************************************************\
 ***  FUNCTION:
 ***     
 *** 
 ***  DESCRIPTION:
 ***
 ***            
 ***  PARAMETERS:
 ***      Type               Name              Range              Description
 ***      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 ***      
 ***  RETURNS:
 ***      
 ***  SETTINGS:
 ***
\******************************************************************************************************************************************************************************************************/
mtxResultInfo mtx_sub_scalar_f64(tMatrix * const pSrc,const float64 scalar)
{
    mtxResultInfo Result = MTX_OPERATION_OK;
    float64 * const pDst= pSrc->val;    
    uint16 eIdx;
    
    for(eIdx=0;eIdx<pSrc->nelem;eIdx++)
    {
        pDst[eIdx] -= scalar;
    }
    
    return Result;
}
/******************************************************************************************************************************************************************************************************\
 ***  FUNCTION:
 ***     
 *** 
 ***  DESCRIPTION:
 ***
 ***            
 ***  PARAMETERS:
 ***      Type               Name              Range              Description
 ***      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 ***      
 ***  RETURNS:
 ***      
 ***  SETTINGS:
 ***
\******************************************************************************************************************************************************************************************************/
mtxResultInfo mtx_add_scalar_f64(tMatrix * const pSrc,const float64 scalar)
{
    mtxResultInfo Result = MTX_OPERATION_OK;
    float64 * const pDst= pSrc->val;    
    uint16 eIdx;
    
    for(eIdx=0;eIdx<pSrc->nelem;eIdx++)
    {
        pDst[eIdx] += scalar;
    }
    
    return Result;
}
/******************************************************************************************************************************************************************************************************\
 ***  FUNCTION:
 ***     
 *** 
 ***  DESCRIPTION:
 ***
 ***            
 ***  PARAMETERS:
 ***      Type               Name              Range              Description
 ***      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 ***      
 ***  RETURNS:
 ***      
 ***  SETTINGS:
 ***
\******************************************************************************************************************************************************************************************************/ 
mtxResultInfo mtx_cpy_f64(tMatrix * const pDst,tMatrix const * const pSrc)
{
    mtxResultInfo Result = MTX_OPERATION_OK;
    float64 * const pDstL = pDst->val;
    float64 const * const pSrcL= pSrc->val;
    uint16 eIdx;

    if(pDst->ncol == pSrc->ncol && pDst->nrow == pSrc->nrow)
    {
        for(eIdx=0;eIdx<pSrc->nelem;eIdx++)
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
/******************************************************************************************************************************************************************************************************\
 ***  FUNCTION:
 ***     
 *** 
 ***  DESCRIPTION:
 ***
 ***            
 ***  PARAMETERS:
 ***      Type               Name              Range              Description
 ***      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 ***      
 ***  RETURNS:
 ***      
 ***  SETTINGS:
 ***
\******************************************************************************************************************************************************************************************************/ 
mtxResultInfo mtx_identity_f64(tMatrix * const pSrc)
{
    mtxResultInfo Result = MTX_OPERATION_OK;
    float64 * const pDst= (float64 *)pSrc->val;    
    const uint8 nCol = pSrc->ncol;
    uint16 eIdx;
    
    if(pSrc->nrow == nCol)
    {   
        pDst[0] = 1;
        
        for(eIdx=1;eIdx < pSrc->nelem;eIdx++)
        {
            const uint16 cmpLeft = (uint16)(eIdx / nCol);
            
            pDst[eIdx] = eIdx < nCol ? 0 : cmpLeft == eIdx % (cmpLeft*nCol) ?  1 : 0;   
        }    
    }
    else
    {
        Result = MTX_SIZE_MISMATCH;
    }

    return Result;
}
/******************************************************************************************************************************************************************************************************\
 ***  FUNCTION:
 ***     
 *** 
 ***  DESCRIPTION:
 ***
 ***            
 ***  PARAMETERS:
 ***      Type               Name              Range              Description
 ***      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 ***      
 ***  RETURNS:
 ***      
 ***  SETTINGS:
 ***
\******************************************************************************************************************************************************************************************************/ 
mtxResultInfo mtx_zeros_f64(tMatrix * const pSrc)
{
    mtxResultInfo Result = MTX_OPERATION_OK;
    float64 * const pDst = (float64 *)pSrc->val;
    uint16 eIdx;    
       
    for(eIdx=0;eIdx<pSrc->nelem;eIdx++)
    {
        pDst[eIdx] = 0;             
    }
    
    return Result;
}
