#include "iostream"
#include "function.h"
#include "mesh.h"
#include <omp.h>
#include "bc.h"
#include "mpi.h"

void apply_BC_init(int mesh_size,int overlapping,double ***K,double **b){
  int my_rank;
  MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
  double start,end,start_total,end_total;
  double **K_local_bc,*b_local_bc;
  K_local_bc = *K;
  b_local_bc = *b;
  //sol = (double *) malloc ((mesh_size+1)*(mesh_size+1)*sizeof(double));
  //start_total = omp_get_wtime();
  //start = omp_get_wtime();
  //#pragma omp parallel for shared(mesh_size,K_global,b_global)
    //set the bottom as D
    for (int i = 0;i< mesh_size+1;i++)
    {
      for(int j = 0;j<(mesh_size/2+overlapping+1)*(mesh_size+1);j++)
      {
        K_local_bc[i][j] = 0.0;
        K_local_bc[j][i] = 0.0;
      }
      K_local_bc[i][i] = 1.0;
      b_local_bc[i] = 0.0;
    }

    //set the top as D
    for(int i = (mesh_size/2+overlapping)*(mesh_size+1);i< (mesh_size/2+1+overlapping)*(mesh_size+1);i++)
    {
      for(int j = 0;j<(mesh_size/2+overlapping+1)*(mesh_size+1);j++)
      {
        K_local_bc[i][j] = 0.0;
        K_local_bc[j][i] = 0.0;
      }
      K_local_bc[i][i] = 1.0;
      b_local_bc[i] = 0.0;
    }

  //end = omp_get_wtime();
  //printf(" - For BC part1: \t %f seconds\n", end - start);

  //start = omp_get_wtime();
//  #pragma omp parallel for shared(mesh_size,K_global,b_global)
  for(int i = mesh_size+1;i<(mesh_size/2+1+overlapping)*(mesh_size+1);i=i+(mesh_size+1))
  {
    for(int j = 0;j<(mesh_size/2+overlapping+1)*(mesh_size+1);j++)
    {
      K_local_bc[i][j] = 0.0;
      K_local_bc[j][i] = 0.0;
    }
    K_local_bc[i][i] = 1.0;
    b_local_bc[i] = 0.0;
  }

  for(int i = mesh_size;i<(mesh_size/2+1+overlapping)*(mesh_size+1);i=i+(mesh_size+1))
  {
    for(int j = 0;j<(mesh_size/2+overlapping+1)*(mesh_size+1);j++)
    {
      K_local_bc[i][j] = 0.0;
      K_local_bc[j][i] = 0.0;
    }
    K_local_bc[i][i] = 1.0;
    b_local_bc[i] = 0.0;
  }



  *K = K_local_bc;
  *b = b_local_bc;
}


void apply_BC_inteface(mesh mesh_in,int mesh_size,int overlapping,double ***K,double **b){
  int my_rank;
  int length_y = mesh_size/2+overlapping+1;
  MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
  double *b_local_inter,**K_local_inter;
  b_local_inter = *b;
  K_local_inter = *K;
  int *dirichlet = new int[mesh_size+1];

  if(my_rank==0){
    for (int i = 0; i < mesh_size+1; i++) {
      dirichlet[i] = i+(length_y-1)*(mesh_size+1);

    }
  }
  else{
    for (int i = 0; i < mesh_size+1; i++) {
      dirichlet[i] = i;
    }
  }

  int tmp_index;
  if(my_rank==0){
    for(int i = (length_y-1)*(mesh_size+1);i<(length_y)*(mesh_size+1);i++){
      for(int j = 0;j<(length_y-1)*(mesh_size+1);j++){
          b_local_inter[j] = b_local_inter[j] - K_local_inter[j][i]*b_local_inter[i];
      }
    }
  }else{
    for(int i = 0;i<(mesh_size+1);i++){
      for(int j = (mesh_size+1);j<(mesh_size+1)*(length_y);j++){
          b_local_inter[j] = b_local_inter[j] - K_local_inter[j][i]*b_local_inter[i];
      }
    }
  }


  for(int i = 0;i<mesh_size+1;i++){
    tmp_index = dirichlet[i];
    for(int j = 0;j<(mesh_size+1)*(length_y);j++){
      K_local_inter[tmp_index][j] = 0;
      K_local_inter[j][tmp_index] = 0;
    }
    K_local_inter[tmp_index][tmp_index] = 1;
  }

  for(int i = mesh_size+1;i<(mesh_size/2+1+overlapping)*(mesh_size+1);i=i+(mesh_size+1))
  {
    for(int j = 0;j<(mesh_size/2+overlapping+1)*(mesh_size+1);j++)
    {
      K_local_inter[i][j] = 0.0;
      K_local_inter[j][i] = 0.0;
    }
    K_local_inter[i][i] = 1.0;
    b_local_inter[i] = 0.0;
  }

  if(my_rank == 0){
    //set bottom as D
    for(int i = 0;i<(mesh_size+1);i++)
    {
      for(int j = 0;j<(mesh_size/2+overlapping+1)*(mesh_size+1);j++)
      {
        K_local_inter[i][j] = 0.0;
        K_local_inter[j][i] = 0.0;
      }
      K_local_inter[i][i] = 1.0;
      b_local_inter[i] = 0.0;
    }
  }
  else
  {
    //set top as 0
    for(int i = (mesh_size/2+overlapping)*(mesh_size+1);i< (mesh_size/2+1+overlapping)*(mesh_size+1);i++)
    {
      for(int j = 0;j<(mesh_size/2+overlapping+1)*(mesh_size+1);j++)
      {
        K_local_inter[i][j] = 0.0;
        K_local_inter[j][i] = 0.0;
      }
      K_local_inter[i][i] = 1.0;
      b_local_inter[i] = 0.0;
    }
  }

  *b = b_local_inter;
  *K = K_local_inter;

}
