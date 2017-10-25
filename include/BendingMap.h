#ifndef BENDINGMAP_H
#define BENDINGMAP_H

#include "ElementMap.h"

class BendingMap:public ElementMap
{
public:
  BendingMap();
  virtual ~BendingMap();  
  BendingMap(double length,double radius,double aperture,double gamma=2000);//mm,m,mm
  void Map();
  void setradius(double);
  void setangle(double);
  void setaperture(double);
  void print();
protected:
  double h,angle,radius;
};

#endif
