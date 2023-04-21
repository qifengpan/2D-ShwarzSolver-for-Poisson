#ifndef point_h
#define point_h
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
    int    global_id;
    int   is_interface_dirichlet=0;
public:

    point(){}//{id = 0;};
    ~ point(){};
    double get_X();
    double get_Y();
    void set_X(double);
    void set_Y(double);

    int get_ID();
    int get_G_ID();
    void set_ID(int,int);
    void set_dirichlet();
    int get_dirichlet();
};

#endif
