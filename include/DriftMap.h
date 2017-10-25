#ifndef DRIFTMAP_H
#define DRIFTMAP_H

#include "ElementMap.h"
class DriftMap:public ElementMap
{
public:
  virtual ~DriftMap();
  DriftMap();
  DriftMap(double length,double gamma);
  const vector<double> getfield(const vector<double> &xyz);
  void Map();
  void print();
};
#endif
