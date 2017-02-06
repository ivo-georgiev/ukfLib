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
int MtxSumDownDiag(int*A ,int m, int n);
void MtxTr(int* A,int m,int n);
/*---------------------------------------------*/
/*         Macros definiton                    */
/*---------------------------------------------*/
#define columns(matrix) sizeof(matrix[0])/sizeof(matrix[0][0])
#define rows(matrix)    sizeof(matrix)/sizeof(matrix[0])
#define rowxcol(matrix) sizeof(matrix)/sizeof(matrix[0][0])

/*---------------------------------------------*/
/*         Structure prototype                 */
/*---------------------------------------------*/
struct matrix
{
    int row;
    int col;
    int* val;
};

struct matrix MtxMul1(struct matrix A, struct matrix B,struct matrix C);

/*---------------------------------------------*/
/*             Sample Matrix                   */
/*---------------------------------------------*/
double symMtx[5][5] =
{{ 3.009964742796415,  -1.009481719676593,   1.161774056277037,  -0.291928007896218,  -0.775215823251770},
{ -1.009481719676593,   3.009964742796415,  -1.009481719676593,   1.161774056277037,  -0.291928007896218},
{  1.161774056277037,  -1.009481719676593,   3.009964742796415,  -1.009481719676593,   1.161774056277037},
{ -0.291928007896218,   1.161774056277037,  -1.009481719676593,   3.009964742796415,  -1.009481719676593},
{ -0.775215823251770,  -0.291928007896218,   1.161774056277037,  -1.009481719676593,   3.009964742796415}};

//1.734924996302842      -0.581858997840148   0.669639355449256      -0.168265491890613  -0.446829589119858
//0                       1.634443284249678  -0.379239855780691       0.650904975441148  -0.337680621986338
//0                       0                   1.554903536627710      -0.418003689501540   0.857240820764834
//0                       0                   0                       1.543776893059448  -0.328117294491480
//0                       0                   0                       0                   1.361527478565284
double Px[2][2]={{0.8,0},{0, 0.3}};

double dblC[3][3] = { {10.5,2.17,3.03},   //size 3x3
                      {0.44,0.59,6.89},
                      {7.56,8.17,9.21}};

double dblD[3][3] = { {1.11,29.3,31.2},   //size 3x3
                      {45.3,5.17,6.11},
                      {7.61,88.0,9.34}};
double dblM[3][3];

double sz[2][3] = { {1.11,29.3,31.2},   //size 3x3
                    {45.3,5.17,6.11}};

int C[3][3] = { {1,2,3},   //size 3x3
                {4,5,6},
                {7,8,9}};

int D[3][3] = { {1,2,3},   //size 3x3
                {4,5,6},
                {7,8,9}};

int F[2][3] = { {3,2,1},  //size 2x3
                {6,5,4}};

int Q[3][2] = { {1,2},    //size 3x2
                {4,5},
                {7,8}};

int R[3][1] = { {1},      //size 3x1
                {2},
                {3}};
int CC[3][3];// product matrix
int CR[3][1];// product matrix
int M[3][3]; // destination matrix
int N[3][3]; // destination matrix


void main(void)
{
    struct matrix a1,a2,a12 ;
    int i=0,j=0;
    int *ptr;
    int *ptr1;
    double *ptrDbl;
//    double upL[5][5];

    printf("num cols = %d \n ",columns(sz));
    printf("num rows = %d \n",rows(sz));
    printf("num rows*cols = %d \n ",rowxcol(sz));

    puts("TEST 1 - Cholesky factorization " );
    //(void)MtxCholUpper(&symMtx[0][0], &upL[0][0],5);
    mtxLib_Cholesky_LL_dp(&Px[0][0],2);
    /*-----------------------------------------*/ 
    /*        TEST 0 Matrix product in Double  */
    /*        M=C*D result is 3x3 matrix       */
    /*-----------------------------------------*/ 
    puts("TEST 1 - Matrix product " );
    ptrDbl = doubleMtxMul(&dblC[0][0],&dblD[0][0],&dblM[0][0],3,3);
    puts("dblM=dblC*dblD");
    for(i=0;i<3;i++)
    {
        putchar('|');
        for(j=0;j<3;j++)
        {
            printf("%f ",dblC[i][j]);       
        }
        putchar('|');
        i==1? putchar('*'): putchar(' ');
        putchar('|');
        for(j=0;j<3;j++)
        {
            printf("%f ",dblD[i][j]);       
        }
        putchar('|');
        i==1? putchar('='): putchar(' ');
        putchar('|');
        for(j=0;j<3;j++)
        {
            printf("%f ",*(ptrDbl+3*i+j)); 
        }
        //putchar('|');
        puts("|");
    }
    putchar('\n');


    /*-----------------------------------------*/ 
    /*        TEST 1 Matrix product            */
    /*        M=C*D result is 3x3 matrix       */
    /*-----------------------------------------*/ 
    puts("TEST 1 - Matrix product " );
    ptr = MtxMul(&C[0][0],&D[0][0],&M[0][0],3,3);
    puts("M=C*D");
    for(i=0;i<3;i++)
    {
        putchar('|');
        for(j=0;j<3;j++)
        {
            printf("%8d",C[i][j]);       
        }
        putchar('|');
        i==1? putchar('*'): putchar(' ');
        putchar('|');
        for(j=0;j<3;j++)
        {
            printf("%8d",D[i][j]);       
        }
        putchar('|');
        i==1? putchar('='): putchar(' ');
        putchar('|');
        for(j=0;j<3;j++)
        {
            printf("%8d",*(ptr+3*i+j)); 
        }
        //putchar('|');
        puts("|");
    }
    putchar('\n');

    /*-----------------------------------------*/
    /*       TEST 2 Matrix transponce          */
    /*-----------------------------------------*/
    puts("TEST 2 - Matrix transponce fnc1 " );
    /*   transp operation test function 1 */
    ptr1 = MtxTrDest(&C[0][0],&N[0][0],3,3);
    for(i=0;i<3;i++)
    {
        putchar('|');
        for(j=0;j<3;j++)
        {
            printf("%8d",C[i][j]);       
        }
        putchar('|');
        i==1? putchar('='): putchar(' ');
        putchar('|');

        for(j=0;j<3;j++)
        {           
            printf(" %8d",*(ptr1+i*3+j)); 
        }
        puts("|");
    }
    putchar('\n');
    puts("TEST 2 Matrix transponce fnc2 " );
    /*   transp operation test function 2 */
    MtxTr(&N[0][0],3,3);
    putchar('\n');

    /*-----------------------------------------*/ 
    /*        TEST 3 diagonal sum              */
    /*-----------------------------------------*/ 
    puts("TEST 3 - diagonal sum of N  " );
    
    
    for(i=0;i<3;i++)
    {
       
        putchar('|');

        for(j=0;j<3;j++)
        {           
            printf(" %8d",*(&N[0][0]+i*3+j)); 
        }       
        puts("|");   
    }
    printf("Sum of diag of N = %8d\n",MtxSumDiag(&N[0][0],3, 3));
    putchar('\n');

    /*-----------------------------------------*/ 
    /*        TEST 4 upper diagonal sum        */
    /*-----------------------------------------*/ 
    puts("TEST 4 - upper diagonal sum " );
    for(i=0;i<3;i++)
    {
       
        putchar('|');

        for(j=0;j<3;j++)
        {           
            printf(" %8d",*(&N[0][0]+i*3+j)); 
        }       
        puts("|");   
    }
    printf("\nUpper diag sum = %8d\n",MtxSumUpDiag(&N[0][0],3, 3));
    putchar('\n');

    /*-----------------------------------------*/ 
    /*        TEST 5  down diagonal sum        */
    /*-----------------------------------------*/ 
    puts("TEST 5 - down diagonal sum " ); 
    for(i=0;i<3;i++)
    {
       
        putchar('|');

        for(j=0;j<3;j++)
        {           
            printf(" %8d",*(&N[0][0]+i*3+j)); 
        }       
        puts("|");   
    }
    printf("\n down diag sum of matrix = %8d\n",MtxSumDownDiag(&N[0][0],3, 3));
    putchar('\n');

    /*------------------------------------------*/ 
    /*        TEST 6 number of columns          */
    /*------------------------------------------*/ 
    puts("TEST 6 - number of columns on matrix" ); 
    printf("\n number of columns of F = %d\n",columns(F));
    putchar('\n');

    /*------------------------------------------*/ 
    /*        TEST 7 number of rows             */
    /*------------------------------------------*/ 
    puts("TEST 7 - number of rows on matrix" ); 
    printf("\nnumber of rows of F = %d\n",rows(F));
    putchar('\n');

    /*------------------------------------------*/ 
    /*        TEST 8 rows*col                   */
    /*------------------------------------------*/ 
    puts("TEST 8 - product rows*col on F matrix" ); 
    printf("\n rows*col = %d\n",rowxcol(F));
    putchar('\n');

    /*------------------------------------------*/ 
    /*    TEST 9 Matrix product by using of     */
    /*    struct M=C*R result is 3x1 matrix     */
    /*------------------------------------------*/ 
    puts("TEST 1 - Matrix product by using stucture " );
    //first matrix
    a1.row = rows(C);
    a1.col = columns(C);
    a1.val = &C[0][0];

    //second matrix
    a2.row = rows(R);
    a2.col = columns(R);
    a2.val = &R[0][0];

    //result matrix
    a12.row = 3;
    a12.col = 1;
    a12.val = &CR[0][0];

    MtxMul1(a1,a2,a12);
    for(i=0;i<a12.row;i++)
    {
        putchar('|');
        for(j=0;j<a1.col;j++)
        {
            //printf("\n %8d\n",M[i][j]); 
            printf(" %8d",*(a1.val+i*a1.col+j)); 
        }    
        putchar('|');
        i==1? putchar('*'): putchar(' ');
        putchar('|');
        for(j=0;j<a2.col;j++)
        {
            //printf("\n %8d\n",M[i][j]); 
            printf(" %8d",*(a2.val+i*a2.col+j)); 
        } 
        putchar('|');
        i==1? putchar('='): putchar(' ');
        putchar('|');
        for(j=0;j<a12.col;j++)
        {
            //printf("\n %8d\n",M[i][j]); 
            printf(" %8d",*(a12.val+i*a12.col+j)); 
        }    
        puts("|");
    }
        putchar('\n');
}

//C<-A*B
int* MtxMul(int* A, int* B, int *C, int m,int n)
{
    int i,j,k;
    int sum;
    for(i=0;i<m;i++)
    {        
        for(j=0;j<n;j++)
        {
            sum = 0;
            for(k=0;k<n;k++)
            {
                /* sum = sum+ A[i][k] * B[k][j];*/
                sum += *(A+i*n+k) * *(B+k*n+j);
            }
            *(C+i*n+j) = sum; //C[i][j]             
        }      
    }
    return C;
}

double* doubleMtxMul(double* A, double* B, double* C, int m,int n)
{
    int i,j,k;
    double sum;
    for(i=0;i<m;i++)
    {        
        for(j=0;j<n;j++)
        {
            sum = 0;
            for(k=0;k<n;k++)
            {
                /* sum = sum+ A[i][k] * B[k][j];*/
                sum += *(A+i*n+k) * *(B+k*n+j);
            }
            *(C+i*n+j) = sum; //C[i][j]             
        }      
    }
    return C;
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
int MtxSumUpDiag(int*A ,int m, int n)
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
int MtxSumDownDiag(int*A ,int m, int n)
{
    int i,j,sum=0;
    for(i=1;i<m;i++)
    {
        for(j=0;j<n;j++)
        {
            if(i > j)
            {
                sum += *(A+i*n+j);
            }
        }
    } 
    return sum;
}

struct matrix MtxMul1(struct matrix A, struct matrix B,struct matrix C)
{
    int i,j,k;
    int sum;
    if(A.col != B.row)
    {
       printf("The rows of A and columns of B must be the same");
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
    return C;
}










