#include "iostream"
#include "mesh.h"
#include "function.h"
#include <omp.h>
#include <cmath>
#include <fstream>
#include "mpi.h"

using namespace std;
//--------------------------------------------------------//
int mesh::visit_node_isdirichlet(int id){
  return point_set[id].get_dirichlet();
}

void mesh::set_node_dirichlet(int id){
  point_set[id].set_dirichlet();
}
element mesh::visit_element(int id)
{
  return mesh_element[id];
}
void mesh::generate_mesh(int n,int delta)
{
  double inc;
  //double start,end;
  //generate_point(n)
  mesh_size = n;
  overlapping = delta;
  //linsp = linspace(0,1,n+1);
  inc = (1.0- 0.0)/double(n);
  mesh_element = new element[n*n + 2*n*overlapping];  // This does not parallelize -> The constructor has to be changed

  //printf(" - new element: \t %f seconds\n", end - start);

  generate_element(mesh_size,inc,overlapping);
}


void mesh::generate_element(int n_elem,double inc,int overlapping)
{
  int n_point;
  int stride,offset_id,gama_judge;
  double offset_y;
  int my_rank;

  MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
  if(my_rank == 0){
    offset_id = 0;
    offset_y = 0;
    gama_judge = n_elem/2+overlapping;  //the top row of rank0 is gama boundary
  }
  else{
    offset_id = (n_elem/2 + 1-overlapping) * (n_elem + 1);
    offset_y = inc * double(n_elem/2-overlapping);
    gama_judge = 0; // the first row of rank1 is gama boundary
  }
  n_point= n_elem + 1;
  stride = n_elem * n_elem / 2 + n_elem * overlapping;

  point *point_set = new point[n_point*(n_elem / 2 + 1 + overlapping)];


  //double start,end;
  //start = omp_get_wtime();
  // #pragma omp parallel for shared(n_point,point_set)
  //#pragma omp parallel for shared(n_point)
  for(int i=0;i< n_elem/2+ 1 + overlapping;i++)
  {
    for(int j=0;j< n_point;j++)
    {
      //printf("size of j %ld\n",sizeof(j));
      point_set[i*n_point+j].set_X(double(j)*inc);
      point_set[i*n_point+j].set_Y(double(i)*inc + offset_y);
      point_set[i*n_point+j].set_ID(i*n_point+j,offset_id);

      //printf("y is %lf\n",point_set[i*n_point+j].get_Y() );

    }
  }

  //end = omp_get_wtime();
  //printf(" - For points: \t\t %f seconds\n", end - start);

  // Parallelization of this loop?
  //start = omp_get_wtime();
  //#pragma omp parallel for shared(n_elem,mesh_element,stride,n_point,point_set)
  //#pragma omp parallel for shared(n_elem)
  for(int i=0; i < n_elem/2+overlapping; i++)
  {
    for (int j = 0; j < n_elem; j++) {
      mesh_element[i*n_elem+j].assign_point( point_set,i*n_point+j);
      mesh_element[i*n_elem+j].assign_point( point_set,i*n_point+j+1);
      mesh_element[i*n_elem+j].assign_point( point_set,i*n_point+j+n_point);

      mesh_element[i*n_elem+j+stride].assign_point( point_set,i*n_point+j+1+n_point);
      mesh_element[i*n_elem+j+stride].assign_point( point_set,i*n_point+j+n_point);
      mesh_element[i*n_elem+j+stride].assign_point( point_set,i*n_point+j+1);
      if((i*n_point+j+n_point+1)%n_point == n_elem){
        mesh_element[i*n_elem+j+stride].set_neumann();
      }
      if(i == gama_judge){
        mesh_element[i*n_elem+j+stride].set_gama();
      }
    }
  }

  //end = omp_get_wtime();
  //printf(" - For elements: \t %f seconds\n", end - start);

  //start = omp_get_wtime();
  // #pragma omp parallel for shared(n_elem,mesh_element)
  //#pragma omp parallel for shared(n_elem)
  for (int i= 0; i < n_elem*n_elem + 2 * n_elem * overlapping; i++) {
    mesh_element[i].generate_middlepoint();
    mesh_element[i].generate_localstiff();
    mesh_element[i].generate_localb(mesh_size);
  }

  //end = omp_get_wtime();
  //printf(" - For middle points: \t %f seconds\n", end - start);

}


double mesh::calculate_integral()
{
  double length;
  double integral=0.0;
  double x,y;
  double scale;
  length = fabs(mesh_element[0].visit_node(1).get_X()-mesh_element[0].visit_node(2).get_X());

  scale = 1.0/6.0*(length*length);

  //double start,end;
  //start = omp_get_wtime();
  //#pragma omp parallel for shared(mesh_size,mesh_element,scale) private(x,y) reduction(+:integral)
  //#pragma omp parallel for shared(mesh_size) private(x,y) reduction(+:integral)
  for (int i= 0;i<mesh_size*mesh_size + 2*mesh_size*overlapping;i++)
  {
    //cout << omp_get_thread_num() << " " << i << " " << 2*mesh_size*mesh_size << endl;
    for (int j= 0;j<3;j++)
    {
      x = mesh_element[i].visit_middlepoint_x(j);
      y = mesh_element[i].visit_middlepoint_y(j);
//      std::cout << scale*fun(x,y) << '\n';
      integral = integral + scale*fun(x,y);
    }
  }
  //end = omp_get_wtime();
  //printf(" - For integral: \t %f seconds\n", end - start);
  return integral;
}

element mesh::get_element(int id){
  return mesh_element[id];
}

int mesh::get_mesh_size(){
  return mesh_size;
}

int mesh::get_overlapping(){
  return overlapping;
}
