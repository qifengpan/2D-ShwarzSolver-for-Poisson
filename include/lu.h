#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int LU_kij(double ***A1,double **b1,int n);

int fw_bw_solve(double **A,double *b,double **x1,int n);
