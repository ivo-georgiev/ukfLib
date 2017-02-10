#include<stdio.h>
#include<math.h>
#include <stdlib.h>
#include "mtxLib.h"

/*---------------------------------------------*/
/*           Function Prototype                */
/*---------------------------------------------*/

double* doubleMtxMul(double* A, double* B, double* C, int m,int n);
int* MtxMul(int* A, int* B, int* C, int n,int m);
int* MtxTrDest(int* A,int* B,int m,int n);
int MtxSumDiag(int*A ,int m, int n);
int MtxSumUpDiag(int*A ,int m, int n);
void MtxTr(int* A,int m,int n);

void show_matrix_obj(sMatrixType A);
void show_matrix(double * A, int n,int m);


/*---------------------------------------------*/
/*             Sample Matrix                   */
/*---------------------------------------------*/
double symMtx[5][5] =
{{ 3.009964742796415,  -1.009481719676593,   1.161774056277037,  -0.291928007896218,  -0.775215823251770},
{ -1.009481719676593,   3.009964742796415,  -1.009481719676593,   1.161774056277037,  -0.291928007896218},
{  1.161774056277037,  -1.009481719676593,   3.009964742796415,  -1.009481719676593,   1.161774056277037},
{ -0.291928007896218,   1.161774056277037,  -1.009481719676593,   3.009964742796415,  -1.009481719676593},
{ -0.775215823251770,  -0.291928007896218,   1.161774056277037,  -1.009481719676593,   3.009964742796415}};

double symMtxChol[5][5] =
{{1.734924996302842,     -0.581858997840148,   0.669639355449256,      -0.168265491890613,  -0.446829589119858},
{0                 ,      1.634443284249678,  -0.379239855780691,       0.650904975441148,  -0.337680621986338},
{0                 ,      0                ,   1.554903536627710,      -0.418003689501540,   0.857240820764834},
{0                 ,      0                ,   0                ,       1.543776893059448,  -0.328117294491480},
{0                 ,      0                ,   0                ,       0                ,   1.361527478565284}};

double Identity_5x5[4][4] =
{{1.0,  0     ,0,     0,     },
{0,     1.0   ,0,     0,     },
{0,     0     ,1.0,   0,     },
{0,     0     ,0,     1.0   }};

double TestMatrix_0_4x4[4][4] =
{{3.0,  5.0, -1.0,  -4},
{ 1.0,  4.0, -0.7,  -3},
{ 0,   -2.0,  0,     1},
{-2.0,  6.0,  0,     0.3}};


double TestMatrix_1_3x3[3][3] = 
{{10.5, 2.17, 3.03},
{ 0.44, 0.59, 6.89},
{ 7.56, 8.17, 9.21}};

double TestMatrix_2_3x3[3][3] = 
{ {1.11, 29.3, 31.2},
{45.3, 5.17, 6.11},
{7.61, 88.0, 9.34}};


double TestMatrix_1_2x3[2][3] = 
{ {1.11, 29.3, 31.2},   //size 3x3
  {45.3, 5.17, 6.11}};


void main(void)
{
    sMatrixType myFactMatrix;
    sMatrixType myTestMatx={0,0,NULL};
    sMatrixType myChol={0,0,NULL};
    sMatrixType Im={0,0,NULL};
    sMatrixType oTestMatrix_0_4x4={0,0,NULL};

    mtx_init_f64(&myTestMatx,&TestMatrix_1_2x3[0][0],NROWS(TestMatrix_1_2x3),NCOL(TestMatrix_1_2x3));
    mtx_init_f64(&myChol,&symMtxChol[0][0],NROWS(symMtxChol),NCOL(symMtxChol));
    mtx_init_f64(&Im,&Identity_5x5[0][0],NROWS(Identity_5x5),NCOL(Identity_5x5));
    mtx_init_f64(&oTestMatrix_0_4x4,&TestMatrix_0_4x4[0][0],NROWS(TestMatrix_0_4x4),NCOL(TestMatrix_0_4x4));

    //show_matrix(&TestMatrix_1_2x3[0][0],2,3);
    //show_matrix_obj(myTestMatx);


    (void)mtx_init_f64(&myFactMatrix,&symMtx[0][0],NROWS(symMtx),NCOL(symMtx));
    //show_matrix(&symMtx[0][0],5,5);

    (void)mtx_chol_f64(myFactMatrix);
    //show_matrix_obj(myFactMatrix);

    show_matrix_obj(myChol);
    mtx_transp_f64(&myChol);
    show_matrix_obj(myChol);


    show_matrix_obj(Im);
    show_matrix_obj(oTestMatrix_0_4x4);
    mtx_inv_f64(&oTestMatrix_0_4x4, &Im);
    show_matrix_obj(Im);

}

void show_matrix_obj(sMatrixType A)
{
    int i,j;
    
    for (i = 0; i < A.nrow; i++) 
    {
        for (j = 0; j < A.ncol; j++)
        {
            printf("%2.5f ", A.val[A.ncol*i + j]);
        }
        printf("\n");
    }
    printf("\n");
}

void show_matrix(double * A, int n,int m)
{
    int i,j;
    
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < m; j++)
        {
            printf("%2.5f ", A[m*i + j]);
        }
        printf("\n");
    }
    printf("\n");
}

int* MtxTrDest(int* A,int* B,int m,int n)
{
    int i,j;
    for(i=0;i<m;i++)
    {
        for(j=0;j<n;j++)
        {
            *(B+i*n+j) = *(A+j*n+i); 
        }
    }
    return B;
}

void MtxTr(int* A,int m,int n)
{
    int i,j,temp;
    for(i=0;i<m;i++)
    {
        putchar('[');
        for(j=0;j<n;j++)
        {
            if(i != j && i<j)
            {
                temp = *(A+i*n+j);
                *(A+i*n+j) = *(A+j*n+i);
                *(A+j*n+i) = temp;
            }
            printf(" %8d",*(A+i*n+j)); 
        }
        puts(" ]");
    }
}





int MtxSumDiag(int*A ,int m, int n)
{
    int i,j,sum=0;
    for(i=0;i<m;i++)
    {
        for(j=0;j<n;j++)
        {
            if(i == j)
            {
                sum += *(A+i*n+j);
            }
        }
    } 
    return sum;
}





