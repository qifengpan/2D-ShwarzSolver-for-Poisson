#include "lu.h"
#include <omp.h>
int LU_kij(double ***A1,double **b1,int n)
{
    double **A;
    A = *A1;

    double *b;
    b=*b1;

    double max;
    int i_max;
    double *pivot_help;
    double pivot_b;
    int i,j,k;
//    #pragma omp parallel for shared(A,b,n) private(i,j,k,max,i_max,pivot_help,pivot_b)
    for (k=0;k<n;k++)
    {
        //Search for Pivot element
        max=fabs(A[k][k]);
        i_max=k;
        for (i=k;i<n;i++)
        {
            if (fabs(A[i][k])>max)
            {
                max=fabs(A[i][k]);
                i_max=i;
            }
        }

        if (i_max!=k)
        {
            //Swap rows in matrix
            pivot_help=A[k];
            A[k]=A[i_max];
            A[i_max]=pivot_help;

            //Swap accordingly in rhs
            pivot_b=b[k];
            b[k]=b[i_max];
            b[i_max]=pivot_b;
        }

        //Check if pivoting worked out
        if (fabs(A[k][k])<1e-10)
        {
            printf("\n Achtung!! Null-Pivot A%d %d!!!\n\n",k,k);
        }

        //ij Loop (elimination loop)
        for (i=k+1;i<n;i++)
        {
            A[i][k]=A[i][k]/A[k][k];

            for (j=k+1;j<n;j++)
            {
                //Changed k+1 to k, because we want the column to become zero
                A[i][j]-=A[i][k]*A[k][j];
            }
        }
    }
  *A1=A;
  *b1=b;
  return 0;
}

int fw_bw_solve(double **A,double *b,double **x1,int n)
{
    int i,j,k;
    double *x;
    x=*x1;

    double *y = (double *) malloc (n*sizeof(double));

    //Forward solve
  //  #pragma omp parallel for shared(A,y,b,n,x) private(i,j)
    for (i=0;i<n;i++)
    {
        y[i]=b[i];
        for (j=0;j<i;j++)
        {
            y[i]-=A[i][j]*y[j];
        }
    }
    //Backward solve
    for (i=n-1;i>=0;i--)
    {
        x[i]=y[i];
        for (j=n-1;j>i;j--)
        {
            x[i]-=A[i][j]*x[j];
        }
        x[i]=x[i]/A[i][i];
    }

    *x1=x;

    free(y);

    return 0;
}
