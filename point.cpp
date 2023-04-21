#include "point.h"
using namespace std;

double point::get_X(){return x;}
double point::get_Y(){return y;}
int point::get_ID(){return id;}

void point ::set_X(double x_in) {
  x = x_in;
}
void point ::set_Y(double y_in) {
  y = y_in;
}
void point::set_ID(int ID_in,int offset_id){
    id = ID_in;
    global_id = ID_in + offset_id;
}

void point::set_dirichlet(){
  is_interface_dirichlet = 1;
}

int point::get_dirichlet(){
  return is_interface_dirichlet;
}
