#include "function.h"
#include "tgmath.h"
#include <fstream>
#include <iostream>
using namespace std;

#define PI 3.1415926535

double fun(double x,double y)
{
  return sin(x+y)*sin(2*x*x+4*y);
}
double funA(double x,double y)
{
  if ((pow(x,2)+pow(y,2))<1)
  {
    return 1;
  }
  else
  {
    return 1;
  }
}

/*double funA(double x,double y)
{
  return 1;
}*/




double funexat(double x,double y){
  return sin(PI*x)*sin(PI*y);
}

double funf(double x,double y)
{
  return 2*PI*PI*sin(PI*x)*sin(PI*y);
  //return 0;
  //return sin(x+y*y);
}
double fung(double x,double y)
{
  return PI*cos(PI*x)*sin(PI*y)+PI*sin(PI*x)*cos(PI*y);
//  return 1;
}

int write_data(double *res,int count)
{
  ofstream myfile;
  myfile.open("output.csv");
  if(myfile.is_open())
  {
  for (int i = 0;i<count ;++i)
    {
      myfile << res[i] << ",";
    }

    myfile.close();
  }
  else std::cout << "Unable to open file" << '\n';
  return 0;
}
