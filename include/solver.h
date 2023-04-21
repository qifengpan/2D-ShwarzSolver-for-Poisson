#ifndef solver_h
#define solver_h

#include "mesh.h"

using namespace std;

class solver{
private:
  double **K_local;
  double *sol;
  double *b_local;
  double *recv_vec0;
  double *recv_vec1;

  double **f_B,**f_I;
  double ***S,**g_gama,***F;

  double *lambda;

  public:
    void solve(mesh);

    void calculate_s(double ***,double ***,double **,double **,double **);
    void calculate_g_gama(double *,double ***,double **,double *,double *);
    void inverse(double **,int );
    double** Multi_M(double **,double **);
    double** Minus_M(double **,double **);
    double** Plus_M(double **,double **);
    double* Multi_M_V(double **,double *);
    double* Minus_V(double *,double *);
    double* Plus_V(double *,double *);

    double* calculate_newlambda(double *,double **,double **,double *,double);

    double norm(double*,double*,int);
};

#endif
