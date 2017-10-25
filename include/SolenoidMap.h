#ifndef SOLENOIDMAP_H
#define SOLENOIDMAP_H

#include "ElementMap.h"
#include <Eigen/LU>
class SolenoidMap:public ElementMap
{
public:

  SolenoidMap()=default;
  ~SolenoidMap()=default;
  SolenoidMap(double length,double field,double gamma=2000);
  void Map();
  void setfield(double field);
  void print();
  const vector<double> getfield(const vector<double> &xyz);
  
};

#endif
