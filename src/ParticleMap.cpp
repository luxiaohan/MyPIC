#include "ParticleMap.h"
ParticleMap::ParticleMap()
{/*
  Coodinate.clear();
  for(int i =0;i<6;i++)
    Coodinate.push_back(0.0);
    */
    //co=0;
}
ParticleMap::~ParticleMap(){}


void ParticleMap::setparticle(double x,double dx,double y,double dy,double z,double dz)
{

  particle_x=x;
  particle_y=y;
  particle_z=z;
  particle_dx=dx;
  particle_dy=dy;
  particle_dz=dz;
  /*
  for(int i=0;i<6;i++)
  {
    if(co*6<1900)
      hist[co*6+i]=p[i];
  }*/
  /*
  double p[]={x,y,z,dx,dy,dz};
  vector<double> pp(p,p+6);
  State.push_back(pp);
  */
  /*
  Coodinate.clear();
  Coodinate.push_back(x);
  Coodinate.push_back(dx);
  Coodinate.push_back(y);
  Coodinate.push_back(dy);
  Coodinate.push_back(z);
  Coodinate.push_back(dz);
  */
}

double ParticleMap::getlocationx()
{
return particle_x;
}
double ParticleMap::getlocationy()
{
return particle_y;
}
double ParticleMap::getlocationz()
{
return particle_z;
}
double ParticleMap::getdirectionx()
{
return particle_dx;
}
double ParticleMap::getdirectiony()
{
return particle_dy;
}
double ParticleMap::getdirectionz()
{
return particle_dz;
}
/*
double ParticleMap::getlocationx()
{
return Coodinate[0];
}
double ParticleMap::getlocationy()
{
return Coodinate[1];
}
double ParticleMap::getlocationz()
{
return Coodinate[2];
}
double ParticleMap::getdirectionx()
{
return Coodinate[3];
}
double ParticleMap::getdirectiony()
{
return Coodinate[4];
}
double ParticleMap::getdirectionz()
{
return Coodinate[5];
}
*/
