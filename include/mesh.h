#ifndef mesh_h
#define mesh_h
//
#include <vector>
#include "element.h"
//#include "function.h"
#include "tgmath.h"
using namespace std;

class mesh {
private:
  /* data */
  int mesh_size;
  int overlapping;
  vector<double> linsp;
  point *point_set;
  element *mesh_element;
//  double *sol;

  /* methods */
//  vector<double> linspace(double, double,int);
  void generate_element(int,double,int);

public:

  void generate_mesh(int,int);
  double calculate_integral();
  element visit_element(int);
  int get_mesh_size();
  int get_overlapping();
  element get_element(int id);
  int visit_node_isdirichlet(int);
  void set_node_dirichlet(int);
  //double **K_global;
  //double *b_global;
};



#endif
