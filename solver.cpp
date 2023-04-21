#include "lu.h"
#include "assembling.h"
#include "bc.h"
#include "function.h"
#include "solver.h"
#include <iostream>
//#include "lapacke.h"
#include "mpi.h"
#include <cstring>
#include "math.h"
void solver::solve(mesh mesh_in)
{
  int state_LU,state_solve,write_state;
  int mesh_size,overlapping;
  int my_rank;

  //double **K_global;
  //double *sol;
  //double *b;
  MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
  mesh_size = mesh_in.get_mesh_size();
  overlapping = mesh_in.get_overlapping();

  //apply inition boundary condition
  generate_stiffnessK(mesh_in, &K_local, &b_local);
  apply_BC_init(mesh_size, overlapping, &K_local, &b_local);
  int length_y = mesh_size/2+overlapping+1;
  sol = (double *) malloc (length_y*(mesh_size+1)*sizeof(double));
  state_LU = LU_kij(&K_local,&b_local,length_y*(mesh_size+1));
  state_solve = fw_bw_solve(K_local,b_local,&sol,length_y*(mesh_size+1));

  MPI_Status status;
// do communicate;
//taglist: S -> 0, g_gama -> 1
 recv_vec0 = new double[mesh_size+1];
 recv_vec1 = new double[mesh_size+1];

  if(my_rank == 0){
    //MPI_Recv(b_local+(length_y-1)*(mesh_size+1),mesh_size+1,MPI_DOUBLE,1,11,MPI_COMM_WORLD,&status);
    MPI_Recv(recv_vec0,mesh_size+1,MPI_DOUBLE,1,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
    MPI_Send(sol+(length_y-2*overlapping-1)*(mesh_size+1),mesh_size+1,MPI_DOUBLE,1,0,MPI_COMM_WORLD);
  }
  if(my_rank == 1)
  {
    MPI_Send(sol+(2*overlapping)*(mesh_size+1),mesh_size+1,MPI_DOUBLE,0,11,MPI_COMM_WORLD);
    MPI_Recv(recv_vec1,mesh_size+1,MPI_DOUBLE,0,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
  }


  int count = 0;
  double *sol_b = new double[length_y*(mesh_size+1)];
  for(int i = 0;i<length_y*(mesh_size+1);i++)
  {
    sol_b[i] = 0;
  }


  double *b_temp;
  double **K_temp;
  double *sol_global = new double[(mesh_size+1)*(mesh_size+1)];
  double *sol_global_l = new double[(mesh_size+1)*(mesh_size+1)];
  double time_acc,time_start,time_end;
  while(1){

    for(int i = 0;i<length_y*(mesh_size+1);i++)
    {
      sol_b[i] = sol[i];
    }

    generate_stiffnessK(mesh_in, &K_temp, &b_temp);

    //apply_BC_init(mesh_size, overlapping, &K_temp, &b_temp);
    if(my_rank == 0){
      for(int i = 0;i<mesh_size+1;i++){
        b_temp[i+(length_y-1)*(mesh_size+1)] = recv_vec0[i];
      }
    }else{
      for(int i = 0;i<mesh_size+1;i++){
        b_temp[i] = recv_vec1[i];
      }
    }

    apply_BC_inteface(mesh_in,mesh_size,overlapping,&K_temp, &b_temp);
    time_start = MPI_Wtime();
    state_LU = LU_kij(&K_temp,&b_temp,length_y*(mesh_size+1));
    state_solve = fw_bw_solve(K_temp,b_temp,&sol,length_y*(mesh_size+1));
    time_end = MPI_Wtime();
    time_acc = time_acc + time_end - time_start;
    if(my_rank == 0){
      //MPI_Recv(b_local+(length_y-1)*(mesh_size+1),mesh_size+1,MPI_DOUBLE,1,11,MPI_COMM_WORLD,&status);
      MPI_Recv(recv_vec0,mesh_size+1,MPI_DOUBLE,1,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
      MPI_Send(sol+(length_y-2*overlapping-1)*(mesh_size+1),mesh_size+1,MPI_DOUBLE,1,0,MPI_COMM_WORLD);
    }
    if(my_rank == 1)
    {
      MPI_Send(sol+(2*overlapping)*(mesh_size+1),mesh_size+1,MPI_DOUBLE,0,11,MPI_COMM_WORLD);
      MPI_Recv(recv_vec1,mesh_size+1,MPI_DOUBLE,0,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
    }
    count = count +1 ;
    if(my_rank==0){printf("norm is %lf\n",norm(sol_b,sol,length_y*(mesh_size+1)) );}


    if(my_rank==1){
      MPI_Send(sol+(2*overlapping+1)*(mesh_size+1),(length_y-2*overlapping-1)*(mesh_size+1),MPI_DOUBLE,0,111,MPI_COMM_WORLD);
    }else{
      MPI_Recv(sol_global+(length_y)*(mesh_size+1),(length_y-2*overlapping-1)*(mesh_size+1),MPI_DOUBLE,1,111,MPI_COMM_WORLD,&status);
    }
    if(my_rank ==0){
      for(int i = 0;i<length_y*(mesh_size+1);i++){
        sol_global[i]= sol[i];
      }
    }
    if(my_rank ==0){
      MPI_Send(sol_global,(mesh_size+1)*(mesh_size+1),MPI_DOUBLE,1,222,MPI_COMM_WORLD);
    }
    else{
      MPI_Recv(sol_global,(mesh_size+1)*(mesh_size+1),MPI_DOUBLE,0,222,MPI_COMM_WORLD,&status);
    }
    if(norm(sol_global,sol_global_l,(mesh_size+1)*(mesh_size+1))<0.00000001){
      break;
    }
    for(int i= 0;i<(mesh_size+1)*(mesh_size+1);i++)
    {
      sol_global_l[i] = sol_global[i];
    }
    //free(K_temp);
    //free(b_temp);
  }

  /*double *sol_global = new double[(mesh_size+1)*(mesh_size+1)];
  if(my_rank==1){
    MPI_Send(sol+(2*overlapping+1)*(mesh_size+1),(length_y-2*overlapping-1)*(mesh_size+1),MPI_DOUBLE,0,111,MPI_COMM_WORLD);
  }else{
    MPI_Recv(sol_global+(length_y)*(mesh_size+1),(length_y-2*overlapping-1)*(mesh_size+1),MPI_DOUBLE,1,111,MPI_COMM_WORLD,&status);
  }
  if(my_rank ==0){
    for(int i = 0;i<length_y*(mesh_size+1);i++){
      sol_global[i]= sol[i];
    }
  }*/
  printf("loop exit, time for sovler: %lf \n",time_acc );
  if(my_rank ==0)
  {
      write_data(sol_global,(mesh_size+1)*(mesh_size+1));
  }
  if(my_rank ==0)
  {
      printf("the total iteration is %d\n",count );
  }




}


double solver::norm(double *A,double *B,int length_B){

    double result=0;

    for(int i = 0;i<length_B;i++)
    {
      result += pow(A[i]-B[i],2);
    }
    return sqrt(result);
}
