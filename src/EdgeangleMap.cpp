#include "EdgeangleMap.h"

EdgeangleMap::~EdgeangleMap(){}
EdgeangleMap::EdgeangleMap(){}
EdgeangleMap::EdgeangleMap(double angle,double radius,double gap)//deg,mm,mm
{
this->angle=angle*M_PI/180;
this->radius=radius;
this->gap=gap;
K1_edge=0.45;
K2_edge=2.8;
Map();
}
void EdgeangleMap::Map()
{
P_edge=K1_edge*gap/abs(radius)*(1+sin(angle)*sin(angle))/cos(angle)*(1-K1_edge*K2_edge*gap*tan(angle)/abs(radius));
P_edge=P_edge*M_PI/180;
//xx
R(0,0)=1;
R(0,1)=0;
R(1,0)=tan(angle)/abs(radius)*1000;
R(1,1)=1;
//yy
R(2,2)=1;
R(2,3)=0;
R(3,2)=-tan(angle-P_edge)/abs(radius)*1000;
R(3,3)=1;
//zz
R(4,4)=1;
R(4,5)=0;
R(5,4)=0;
R(5,5)=1;
//
}

void EdgeangleMap::print(){
cout<<"The transfer matrix of edga angle(edge angle = "<<angle/M_PI*180<<" deg,curvature radius = "<<radius<<" mm,pole gap = "<<gap<<" mm) is\n"<<R<<endl;
//cout<<angle<<"   "<<P<<endl;
}

