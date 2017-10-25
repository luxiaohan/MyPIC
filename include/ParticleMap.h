#ifndef PARTICLEMAP_H
#define PARTICLEMAP_H

#include <vector>
using namespace std;
class ParticleMap
{
public:
  virtual ~ParticleMap();
  ParticleMap();
  void setparticle(double x,double dx,double y,double dy,double z,double dz);
  double getlocationx();
  double getlocationy();
  double getlocationz();
  double getdirectionx();
  double getdirectiony();
  double getdirectionz();

  int lost=-1;

protected:
  //double hist[1300];
  //int co=0;
  //vector<vector<double>  >	 State;
  //vector<double>			 Coodinate;
  double particle_x,particle_dx,particle_y,particle_dy,particle_z,particle_dz;
  double particle_beta,particle_gamma,particle_energy,particle_static_energy;
};
#endif
