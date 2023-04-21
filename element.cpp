#include "element.h"
#include "function.h"
using namespace std;

void element::assign_point(point *point_in,int pos) {
  /* code */

  point_tri.push_back(&point_in[pos]);
}

void element::set_neumann()
{
  is_neumann = 1;
}
void element::set_gama(){
  is_gama = 1;
}

point element::visit_node(int id)
{
  return *point_tri[id];
}

double element::visit_middlepoint_x(int id)
{
  return m_point_x[id];
}
double element::visit_middlepoint_y(int id)
{
  return m_point_y[id];
}
void element::generate_middlepoint()
{
//  std::vector<double> m_point_x(3,0);
//  std::vector<double> m_point_y(3,0);
  m_point_x.push_back((point_tri[0]->get_X()+point_tri[1]->get_X())/2.0);
  m_point_y.push_back((point_tri[0]->get_Y()+point_tri[1]->get_Y())/2.0);
  m_point_x.push_back((point_tri[1]->get_X()+point_tri[2]->get_X())/2.0);
  m_point_y.push_back((point_tri[1]->get_Y()+point_tri[2]->get_Y())/2.0);
  m_point_x.push_back((point_tri[0]->get_X()+point_tri[2]->get_X())/2.0);
  m_point_y.push_back((point_tri[0]->get_Y()+point_tri[2]->get_Y())/2.0);
}

/*void element::generate_xi()
{
  xi1 = point_tri[1]->get_X()- point_tri[2]->get_X();
  xi2 = point_tri[2]->get_X()- point_tri[0]->get_X();
  eta1 = point_tri[1]->get_Y()- point_tri[2]->get_Y();
  eta2 = point_tri[2]->get_Y()- point_tri[0]->get_Y();
}*/
void element::generate_localstiff()
{
  double integral_local=0.0;

//  K_local = {{2,-1,-1},{-1,1,0},{-1,0,1}};
  //calculate integral of A(x),since integral requires length^2,while derivative has 1/length^2, it is cancelled

  for (int i=0;i<3;i++)
  {
    integral_local = integral_local + 1.0/6.0*funA(m_point_x[i],m_point_y[i]);
  }
  for (int i=0;i<3;i++)
    {
      for(int j=0;j<3;j++)
      {
        K_local[i][j] = K_local[i][j] * integral_local;
      }
    }
}

void element::generate_localb(int mesh_size)
{
  double det,len_e,mean_y,mean_x;
  double integral_local=0.0;
  det = pow(point_tri[1]->get_X() - point_tri[0]->get_X(),2);
  len_e = fabs(point_tri[1]->get_X() - point_tri[0]->get_X());

/*  for (int i = 0;i<3;i++)
  {
    if (point_tri[i]->get_ID()%(mesh_size+1)==mesh_size)
    {
      integral_local += 1.0/6.0* fung(m_point_x[i],m_point_y[i]);
    }
    else
    {
      integral_local += 1.0/6.0* funf(m_point_x[i],m_point_y[i]);
    }
  }*/
  for (int i = 0;i<3;i++)
  {
    b[i] = det  * funf(point_tri[i]->get_X(),point_tri[i]->get_Y()) *1.0/6.0;
  }

  if (is_neumann == 1){
    mean_y = (point_tri[0]->get_Y() + point_tri[2]->get_Y())/2;
    mean_x = (point_tri[0]->get_X() + point_tri[2]->get_X())/2;
    b[0] = b[0] + 0.5*len_e*fung(mean_x,mean_y);
    b[2] = b[2] + 0.5*len_e*fung(mean_x,mean_y);
  }

}

double element :: get_localmatrix(int id1,int id2)
{
  return K_local[id1][id2];
}

double element :: get_localb(int id)
{
  return b[id];
}
