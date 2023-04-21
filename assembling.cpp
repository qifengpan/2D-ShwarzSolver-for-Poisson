#include "assembling.h"
#include "function.h"
#include <omp.h>
#include <iostream>
#include "mpi.h"
void generate_stiffnessK(mesh mesh_in, double ***K, double **b)
{
  int id1,id2;
  double start,end;
  int mesh_size = mesh_in.get_mesh_size();
  int overlapping = mesh_in.get_overlapping();
  double **K_local_as,*b_local_as;

  int length_y = mesh_size/2+overlapping+1;
  K_local_as = (double**)malloc(sizeof(double*)*length_y*(mesh_size+1));
  b_local_as = (double*)malloc(sizeof(double*)*length_y*(mesh_size+1));
  for (int i=0;i<length_y*(mesh_size+1);i++)
  {
    K_local_as[i] = (double*)malloc(sizeof(double)*length_y*(mesh_size+1));
  }
  //assemnbling the global matrix
  //start = omp_get_wtime();
  //#pragma omp parallel for shared(mesh_size) private(id1,id2)
  for(int i=0;i<2*(length_y-1)*mesh_size;i++)
  {
    for (int j = 0;j<3;j++)
    {
      for (int k = 0;k<3;k++)
      {
        id1 = mesh_in.get_element(i).visit_node(j).get_ID();
        id2 = mesh_in.get_element(i).visit_node(k).get_ID();

        K_local_as[id1][id2] += mesh_in.get_element(i).get_localmatrix(j,k);

      }
      b_local_as[id1]      += mesh_in.get_element(i).get_localb(j);
    }
  }

  *K = K_local_as;
  *b = b_local_as;

  //end = omp_get_wtime();
  //printf("Assemnbling: \t \t %f seconds\n", end - start);
  /*for (int i = 0;i<(mesh_size+1)*(mesh_size+1);i++)
  {
    for (int j= 0;j<(mesh_size+1)*(mesh_size+1);j++)
    {
      printf("%f    ",K_global[i][j] );
    }
    printf("\n" );
  }*/
}
