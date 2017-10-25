#include "RFQMap.h"

RFQMap::~RFQMap()
{
  clear();
}

RFQMap::RFQMap(){}

RFQMap::RFQMap(double length,double voltage,double phase,double radius,double gamma,double frequency,double AAcc,double AFocus,char particletype)
{
  initial(length,voltage,phase,radius,gamma,frequency,AAcc,AFocus,particletype);
}


void RFQMap::initial(double length,double voltage,double phase,double radius,double gamma,double frequency,double AAcc,double AFocus,char particletype)
{
  N=int(length/step);

  dR.resize(N);
  step=length/double(N);
  this->voltage=voltage;
  setlength(length);
  setgamma(gamma);
  this->phase=phase*M_PI/180;
  this->AAcc=AAcc;
  this->AFocus=AFocus;
  omiga=2*M_PI*frequency;
  this->radius=radius;
  this->particletype=particletype;
}
void RFQMap::clear()
{
}

void RFQMap::settype(int pretype,int type,int nexttype)
{
  this->type=type;
  this->pretype=pretype;
  this->nexttype=nexttype;
}


void RFQMap::Map()
{
  gammain=gamma;
  gammaout=gammain;
  tsyn=step/(2*beta*C_light);
  zsyn=step/2;
  R=MatrixXd::Identity(6,6);

  for(int i=0;i<N;i++)
  {
    dR[i]=MatrixXd::Identity(6,6);
    gammain=gammaout;
    if(particletype=='e')
    {
      Q=_Q;
      Mass=_MASS;
    }
    gammaout=gammain+abs(Q/(Mass*C_light*C_light))*M_PI*AAcc*voltage/length/2*sin(omiga*tsyn+phase)*sin(M_PI*zsyn/length)*step;
    betaout=sqrt(1-1/gammaout/gammaout);
    gammasyn=(gammain+gammaout)/2;
    betasyn=sqrt(1-1/gammasyn/gammasyn);
    if(type==2||type==-2)
    {
      C1=1;
      C2=sin(M_PI*zsyn/length);
      C3=sin(M_PI*zsyn/length);
      S=-sign(type);
    }
    else if(type==3)
    {
      C1=1/4*(3*cos(M_PI*zsyn/length/2-M_PI/2)+cos(3*(M_PI*zsyn/length/2-M_PI/2)));
      C2=0;
      C3=sin(M_PI*zsyn/length);
      S=-sign(nexttype);
    }
    else if(type==-3)
    {
      C1=1/4*(3*cos(M_PI*zsyn/length/2)+cos(3*M_PI*zsyn/length/2));
      C2=0;
      C3=sin(M_PI*zsyn/length);
      S=-sign(pretype);
    }
    else if(type==4)
    {
      C1=1;
      C2=1/2*cos(M_PI*zsyn/length)+1/2;
      C3=sin(M_PI*zsyn/length)/2;
      S=-sign(nexttype);
    }
    else if(type==-4)
    {
      C1=1;
      C2=-1/2*cos(M_PI*zsyn/length)+1/2;
      C3=sin(M_PI*zsyn/length)/2;
      S=-sign(pretype);
    }
    kx1=-abs(Q/(2*Mass*C_light*C_light))*step/gammasyn/betasyn/betasyn*cos(omiga*tsyn+phase)*(S*voltage*AFocus*C1/radius/radius-pow((M_PI/length),2)*AAcc*voltage*C2/4);

    kx2=1-abs(Q/(Mass*C_light*C_light))*step*AAcc*voltage/gammasyn/betasyn/betasyn/2*pow((M_PI/length),2)*C3*cos(omiga*tsyn+phase);
    /*cout<<C3<<endl<<cos(omiga*tsyn+phase)<<endl<<step*AAcc*voltage/gammasyn/betasyn/betasyn/2*pow((M_PI/length),2)<<endl;
      cout<<"mark"<<endl;
      cout<<AAcc<<endl;*/
    ky1=-abs(Q/(2*Mass*C_light*C_light))*step/gammasyn/betasyn/betasyn*cos(omiga*tsyn+phase)*(-S*voltage*AFocus*C1/radius/radius-pow((M_PI/length),2)*AAcc*voltage*C2/4);

    ky2=kx2;

    kz1=-abs(Q/(Mass*C_light*C_light))*step*AAcc*voltage/gammasyn/betasyn/betasyn/2*(M_PI/length)*C3*cos(omiga*tsyn+phase);

    kz2=ky2;
    cout<<"k   "<<i<<endl<<kx1<<endl<<kx2<<endl<<ky1<<endl<<ky2<<endl<<kz1<<endl<<kz2<<endl;
    Rt1<<1,step/2/1000,0,1;
    Rt2<<1,0,kx1*1000,kx2;
    Rt3<<1,step/2/1000,0,1;
    // cout<<"Rt2x"<<endl<<Rt2<<endl;
    dR[i].block<2,2>(0,0)=Rt1*Rt2*Rt3;

    Rt1<<1,step/2/1000,0,1;
    Rt2<<1,0,ky1*1000,ky2;
    Rt3<<1,step/2/1000,0,1;
    //   cout<<"Rt2y"<<endl<<Rt2<<endl;
    dR[i].block<2,2>(2,2)=Rt1*Rt2*Rt3;

    Rt1<<1,step/2/gammain/gammain/1000,0,1;
    Rt2<<1,0,kz1*1000,kz2;
    Rt3<<1,step/2/gammaout/gammaout/1000,0,1;
    //    cout<<"Rtz"<<endl<<Rt2<<endl;
    dR[i].block<2,2>(4,4)=Rt1*Rt2*Rt3;

    tsyn=tsyn+step/betaout/C_light;
    zsyn=zsyn+step;
    //cout<<i<<endl<<dR[i]<<endl;

  }
  for(int i=0;i<N;i++)
  {
    R=dR[i]*R;
  }
}

void RFQMap::print()
{
  cout<<"the transfer matrix of RFQ is\n"<<R<<endl;
}


