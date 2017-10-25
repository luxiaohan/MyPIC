#include "ThinlensMap.h"
ThinlensMap::ThinlensMap(double focallengthx,double focallengthy)//mm
{
  fx=focallengthx;
  fy=focallengthy;
  Map();
}
void ThinlensMap::Map()
{
//xx
R(0,0)=1;
R(0,1)=0;
R(1,0)=-1/fx*1000;
R(1,1)=1;
//yy
R(2,2)=1;
R(2,3)=0;
R(3,2)=-1/fy*1000;
R(3,3)=1;
//zz
R(4,4)=1;
R(4,5)=0;
R(5,4)=0;
R(5,5)=1;
//zx
//xz
}

//print
void ThinlensMap::print()
{
cout<<"The transfer matrix of Thin lens(fx = "<<fx<<" mm, fy = "<<fy<<" mm) is\n"<<R<<endl;
cout<<"The twiss matrix of Thin lens(fx = "<<fx<<" mm, fy = "<<fy<<" mm) is\n"<<R2<<endl;
}

