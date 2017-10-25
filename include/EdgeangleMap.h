#ifndef EDGEANGLEMAP_H
#define EDGEANGLEMAP_H

#include "ElementMap.h"
class EdgeangleMap:public ElementMap
{
public:
  virtual ~EdgeangleMap();
  EdgeangleMap();
  EdgeangleMap(double angle,double radiu,double gap);//deg,mm,mm
  void Map();
  void print();
protected:
  double angle,radius,P_edge,gap,K1_edge,K2_edge;
};

#endif
