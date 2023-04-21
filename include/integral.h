#ifndef integral_h
#define integral_h
//
#include<vector>
#include "tgmath.h"
using namespace std;

class point
{
private:
    double x;
    double y;
    int    id;
public:

    point(){}//{id = 0;};
    ~ point(){};
    double get_X();
    double get_Y();
    void set_X(double);
    void set_Y(double);

    int get_ID();
    void set_ID(int);

};

class element {
private:
  /* data */
  std::vector<point*> point_tri;
  std::vector<double> m_point_x;
  std::vector<double> m_point_y;
  double K_local[3][3] = {{2.0,-1.0,-1.0},{-1.0,1.0,0.0},{-1.0,0.0,1.0}};
  double b[3];
//  double length;
  /* methods*/

public:
  element () {}//{point_tri = new point[3];m_point = new point[3];};

  virtual ~element (){};
  void set_point(point *);
  void assign_point(point *,int);
  point visit_node(int);
  double visit_middlepoint_x(int);
  double visit_middlepoint_y(int);
  void generate_middlepoint();

  void generate_localstiff();
  void generate_localb(int);
  double get_localmatrix(int,int);
  double get_localb(int);
};

class mesh {
private:
  /* data */
  int mesh_size;
  vector<double> linsp;
  element *mesh_element;
  double **K_global;
  double *b_global;
  double *sol;
//  point *ponit_set;
  /* methods */
  vector<double> linspace(double, double,int);
  void generate_element(vector<double> v);

public:
  void generate_stiffnessK();
  void generate_mesh(int);
  double calculate_integral();
  element visit_element(int);
  void apply_BC();
  void solve();

};

double fun(double x,double y);
double funA(double,double);
double funf(double,double);
double fung(double,double);
int write_data(double *,int);


#endif
