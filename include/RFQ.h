#ifndef RFQ_H
#define RFQ_H
#include "ElementMap.h"
#include <string>

typedef vector<double> vd;
class RFQ:public ElementMap
{
  public:
    virtual ~RFQ();
    RFQ(const string &path,double position,double length,double frequency);
    RFQ& operator=(const RFQ &rh);
    virtual const vector<double> getfield(const vector<double> &xyz,double time);
    double getradius(double z);

    MatrixXd getmap() const {return MatrixXd::Identity(6,6);}
    MatrixXd getmap()       {return MatrixXd::Identity(6,6);}
  protected:
    vector<vector<double> > factor; 
    vector<double> cell_length;
    vd cell_position;
    vd cell_position_end;
    vd radius;
    vd UL;
    vd mod;
    vd aperture;
  private:
    double ERR=1.0e-16;
    int getEleNum(double z);
    double position,frequency;
    //RFQ::getField private variable
    double syn_phase=0*M_PI/2*2;
};





#endif
