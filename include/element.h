#ifndef element_h
#define element_h
//
#include<vector>
#include "point.h"
//#include "function.h"
#include "tgmath.h"
using namespace std;

class element {
private:
  /* data */
  std::vector<point*> point_tri;
  std::vector<double> m_point_x;
  std::vector<double> m_point_y;
  double K_local[3][3] = {{2.0,-1.0,-1.0},{-1.0,1.0,0.0},{-1.0,0.0,1.0}};
  double b[3];
  int is_neumann;
  int is_gama;
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

  void set_neumann();
  void set_gama();

};


#endif
