#include <iostream>
#include "mesh.h"
#include "solver.h"
//#include <omp.h>
#include "mpi.h"
using namespace std;

int main(int argc, char *argv[]) {
  /* code */
  int  i, size, my_rank;
  MPI_Init(&argc , &argv );
  MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
  MPI_Comm_size(MPI_COMM_WORLD,&size);

  int n;
  double integral;
  if (atoi(argv[1])>1)
    {
      n = atoi(argv[1]);
    }
  else
    {
      cout<<"input number should be positiv"<<endl;
    }

  int delta = atoi(argv[2]); // the overlapping grad.

  //double TOTALstart,TOTALend;
  //TOTALstart = omp_get_wtime();

  mesh trimesh;

  //double start,end;
  //start = omp_get_wtime();
  trimesh.generate_mesh(n,delta);
  //end = omp_get_wtime();
  //printf("Generate mesh: \t\t %f seconds\n", end - start);

//  std::cout << trimesh.visit_element(3).visit_middlepoint(2).get_X() << '\n';
  //start = omp_get_wtime();
  //integral = trimesh.calculate_integral();
  //end = omp_get_wtime();
  //printf("Calculate integral: \t %f seconds\n", end - start);

  /*trimesh.generate_stiffnessK();
  trimesh.apply_BC();
  start = omp_get_wtime();
  trimesh.solve();
  end = omp_get_wtime();
  printf("Solving FEM: \t %f seconds\n", end - start);
  TOTALend = omp_get_wtime();
  printf("TOTAL: \t\t\t %f seconds\n\n", TOTALend - TOTALstart);
  printf("Final result of integral: \t %f \n", integral);
  printf("Successfully excuting FEM, check result in output.dat\n");*/
  double tstart = MPI_Wtime();
  solver solver_tri;
  solver_tri.solve(trimesh);
  double tend = MPI_Wtime();
  printf("time consume by solver(assembling included): %lf\n",tend-tstart );
  MPI_Finalize();
  return 0;
}
