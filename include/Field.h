#ifndef FIELD_H
#define FIELD_H

#include <vector>
#include <string>
#include <ElementMap.h>

using std::vector;
using std::string;
using v3d= vector<vector<vector<double> > > ;
class Field:public ElementMap
{
public:
  vector<double> getfield(vector<double> &xyz);
  double getfield(int);
  void read(const string & );
  void addField(double);
  
protected:
  //1D field
  vector<double> 	 field;
  //3D field
  v3d Ex;
  v3d Ey;
  v3d Ez;
  v3d Bx;
  v3d By;
  v3d Bz;
  
};


#endif
