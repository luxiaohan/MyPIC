#include "BeamMap.h"
#include <fstream>
#include <iostream>
BeamMap::~BeamMap()
{
  clear();
}
BeamMap::BeamMap(){}
BeamMap::BeamMap(char particletype,double gamma,int number)
{
  initial(particletype,gamma,number);
}

void BeamMap::initial(char particletype,double gamma,int number)
{
  if(particletype=='e')
  {
    beam_staticenergy=_MASS*C_light*C_light;
  }
  else if(particletype=='p')
  {
    beam_staticenergy=_MASS*C_light*C_light;
  }
  else
  {
    cout<<"error! particletype not defined"<<endl;
  }
  beam_energy=beam_staticenergy*gamma;
  setgamma_beam_energy(gamma);
  particleLIVEnumber=number;
  particlenumber=number;
  clear();
  particle	.resize(number);
  eGamma	.resize(number);
  eBeta   	.resize(number);
  eBetaGamma	.resize(number);
  phi		.resize(number);
  dw		.resize(number);
}
void BeamMap::clear()
{
  particle.clear();
}

void BeamMap::settwissx(double alpha,double beta)
{
  if(beta>=0)
  {
    TAlphax=alpha;
    TBetax=beta;
    TGammax = (1+TAlphax*TAlphax)/beta;
  }
  else{cout<<"Error! betax= "<<beta<<endl;}
}
void BeamMap::settwissy(double alpha,double beta)
{
  if(beta>=0)
  {
    TAlphay=alpha;
    TBetay=beta;
    TGammay = (1+TAlphay*TAlphay)/beta;
  }
  else{cout<<"Error! betay= "<<beta<<endl;}
}

void BeamMap::setbeta_beam_energy(double a)
{
  if(a>0&&a<1)
  {
    beta_beam_energy=a;
    gamma_beam_energy=1/sqrt(1-beta_beam_energy*beta_beam_energy);
    beam_energy = gamma_beam_energy * beam_staticenergy;
  }
  else if(beta_beam_energy=1)
  {
    beta_beam_energy=a;}
  else
  {
    cout<<"beta_energy wrong!"<<endl;
  }
}

void BeamMap::setgamma_beam_energy(double a)
{
  if(a>=1)
  {
    gamma_beam_energy=a;
    beta_beam_energy=sqrt(1-1/(gamma_beam_energy*gamma_beam_energy));
    beam_energy= gamma_beam_energy* beam_staticenergy;
  }
  else
  {
    cout<<"gamma_energy wrong!"<<endl;
  }
}
void BeamMap::caculate_emittance(const vector<vector<double> > &xyz,const vector<int> &particleLost,double frq)
{
  vector<double> total(6,0.0),sigma(6,0),sigmaxdx(3,0),emit(3);
  particleLIVEnumber	=0;
  eGammaTotal		=0;
  eBetaGammaTotal	=0;
  phiTotal		=0;
  eGammaSigma		=0;
  eBetaGammaSigma	=0;
  phiSigma		=0;
  PhiGamma		=0;
  PhiBetaGamma		=0;
  for(int i=0;i<particlenumber;i++)
  {
    if(particleLost[i]) continue;
    eBeta [i]    = xyz[i][5]/C_light;
    eGamma[i]	 = 1.0/sqrt(  1-pow(eBeta[i],2)  );
    eBetaGamma[i]= eBeta[i]*eGamma[i];
  }
  
  for(int i=0;i<particlenumber;i++)
  {
    if(particleLost[i]) continue;
    for(int j=0;j<6;++j)
    {
      total[j]+=xyz[i][j];
    }
    eGammaTotal		+=eGamma[i];
    eBetaGammaTotal	+=eBetaGamma[i];

    ++particleLIVEnumber;
  }
  eBetaTotal	= total[5]/C_light;

  if(particleLIVEnumber==0) 
  {
    cout<<"No particle LIVE, send in caculate_emittance"<<endl;
    return;
  }
  for(int i=0;i<6;++i)
  {
    particle_ave[i]	=	total[i]/particleLIVEnumber;
  }
  eBetaAverage		= 	eBetaTotal	/particleLIVEnumber;
  eGammaAverage		=	eGammaTotal	/particleLIVEnumber;
  eBetaGammaAverage	=	eBetaGammaTotal	/particleLIVEnumber;
  dwwmax=0;
  dppmax=0;
  dbbmax=0;
  for(int i=0;i<particlenumber;i++)
  {
    if(particleLost[i]) continue;
    for(int j=0;j<6;j=j+2)
    {
      sigma[j]		+=pow( xyz[i][j]-particle_ave[j]    	,2);
      sigma[j+1]	+=pow((xyz[i][j+1]-particle_ave[j+1])/xyz[i][5],2);
    }
    for(int j=0;j<3;++j)
    {
      sigmaxdx[j]	+=(xyz[i][2*j]-particle_ave[2*j])	*	(xyz[i][2*j+1]-particle_ave[2*j+1])	/	xyz[i][5];
    }
    eGammaSigma		+=pow( eGamma[i]-eGammaAverage	,2);
    PhiGamma		+=  (eGamma[i]-eGammaAverage)	*	(xyz[i][4]  -particle_ave[4]);
    eBetaGammaSigma	+=pow( eBetaGamma[i]-eBetaGammaAverage	,2);
    PhiBetaGamma	+=  (eBetaGamma[i]-eBetaGammaAverage)*	(xyz[i][4]  -particle_ave[4]);	
    if(abs(eGamma[i]-eGammaAverage)/(eGammaAverage-1)>dwwmax) 
    {
      dwwmax=abs(eGamma[i]-eGammaAverage)/(eGammaAverage-1);
    }
    if(abs(eBetaGamma[i]-eBetaGammaAverage)/(eBetaGammaAverage)>dppmax) 
    {
      dppmax=abs(eBetaGamma[i]-eBetaGammaAverage)/(eBetaGammaAverage);
    }
    if(abs(eBeta[i]-eBetaAverage)/(eBetaAverage)>dbbmax) 
    {
      dbbmax=abs(eBeta[i]-eBetaAverage)/(eBetaAverage);
    }
  }
  sigmaz	=	eGammaSigma	/ pow(eGammaAverage-1	,2)/particleLIVEnumber;
  eGammaSigma	=	eGammaSigma	/ pow(eGammaAverage	,2);
  eBetaGammaSigma=	eBetaGammaSigma	/ pow(eBetaGammaAverage	,2);
  PhiGamma	=	PhiGamma   	/    (eGammaAverage	)	*   (1/particle_ave[5]*frq*360  );
  PhiBetaGamma	=	PhiBetaGamma	/    (eBetaGammaAverage)	*   (1/particle_ave[5]*frq*360  );
  phiSigma	=	sigma[4]   					*pow(1/particle_ave[5]*frq*360,2);

  

  /*
  if(Ebeta>0&&Ebeta<=1)
  {
    for(int j=1;j<6;j=j+2)
    {
      sigma[j]		/=Ebeta*C_light*Ebeta*C_light;
    }
    for(int j=0;j<3;++j)
    {
      sigmaxdx[j]	/=Ebeta*C_light;
    }
  }
  */
  for(int j=0;j<6;++j)
  {
    sigma[j]	/=	particleLIVEnumber;
  }
  for(int j=0;j<3;++j)
  {
    sigmaxdx[j]	/=	particleLIVEnumber;
  }
  phiSigma	 /=particleLIVEnumber;
  eGammaSigma	 /=particleLIVEnumber;
  PhiGamma	 /=particleLIVEnumber;
  eBetaGammaSigma/=particleLIVEnumber;
  PhiBetaGamma	 /=particleLIVEnumber;

  for(int i=0;i<3;++i)
  {
    emit[i]=sqrt(sigma[2*i]*sigma[2*i+1]-sigmaxdx[i]*sigmaxdx[i]);
  }
  emittancex	=emit[0];
  emittancey	=emit[1];
  x_sigma	=sigma[0];
  dx_sigma	=sigma[1];
  y_sigma	=sigma[2];
  dy_sigma	=sigma[3];
  z_sigma	=phiSigma;
  dz_sigma	=eBetaGammaSigma;
  //dz_sigma	=eGammaSigma;
  //sigmaz	=eGammaSigma;
  sigmadz	=sigma[5];
  xdx_sigma	=sigmaxdx[0];
  ydy_sigma	=sigmaxdx[1];
  zdz_sigma	=PhiBetaGamma;
  //zdz_sigma	=PhiGamma;
  emittancez	=sqrt(z_sigma*dz_sigma-zdz_sigma*zdz_sigma);//*(eGammaAverage-1);
  TBetax	=x_sigma/emittancex;
  TBetay	=y_sigma/emittancey;
  TBetaz	=z_sigma/emittancez;
  TGammax	=dx_sigma/emittancex;
  TGammay	=dy_sigma/emittancey;
  TGammaz	=dz_sigma/emittancez;
  if(TBetax*TGammax>1)
    TAlphax=(xdx_sigma<0?1:-1)*sqrt(TBetax*TGammax-1);
  else
    TAlphax=0;
  if(TBetay*TGammay>1)
    TAlphay=(ydy_sigma<0?1:-1)*sqrt(TBetay*TGammay-1);
  else
    TAlphay=0;
  if(TBetaz*TGammaz>1)
    TAlphaz=(zdz_sigma<0?1:-1)*sqrt(TBetaz*TGammaz-1);
  else
    TAlphaz=0;

  dz_sigma	=eGammaSigma;
  zdz_sigma	=PhiGamma;
  emittancez	=sqrt(z_sigma*dz_sigma-zdz_sigma*zdz_sigma);//*(eGammaAverage-1);
}

void BeamMap::caculate_emittance()
{
  double x_total=0,dx_total=0;
  x_sigma =0;
  dx_sigma=0;
  xdx_sigma=0;
  double y_total=0,dy_total=0;
  y_sigma=0;
  dy_sigma=0;
  ydy_sigma=0;
  double z_total=0,dz_total=0;
  z_sigma=0;
  dz_sigma=0;
  zdz_sigma=0;
  double x_average,y_average,z_average,dx_average,dy_average,dz_average;
  for(int i=0;i<particleLIVEnumber;i++)
  {
    x_total	+=particle[i].getlocationx();
    dx_total	+=particle[i].getdirectionx();
    y_total	+=particle[i].getlocationy();
    dy_total	+=particle[i].getdirectiony();
    z_total	+=particle[i].getlocationz();
    dz_total	+=particle[i].getdirectionz();
  }
  x_average=x_total/particleLIVEnumber;
  dx_average=dx_total/particleLIVEnumber;
  y_average=y_total/particleLIVEnumber;
  dy_average=dy_total/particleLIVEnumber;
  z_average=z_total/particleLIVEnumber;
  dz_average=dz_total/particleLIVEnumber;
  for(int i=0;i<particleLIVEnumber;i++)
  {
    x_sigma=x_sigma+(particle[i].getlocationx()-x_average)*(particle[i].getlocationx()-x_average);
    dx_sigma=dx_sigma+(particle[i].getdirectionx()-dx_average)*(particle[i].getdirectionx()-dx_average);
    xdx_sigma=xdx_sigma+(particle[i].getlocationx()-x_average)*(particle[i].getdirectionx()-dx_average);
    y_sigma=y_sigma+(particle[i].getlocationy()-y_average)*(particle[i].getlocationy()-y_average);
    dy_sigma=dy_sigma+(particle[i].getdirectiony()-dy_average)*(particle[i].getdirectiony()-dy_average);
    ydy_sigma=ydy_sigma+(particle[i].getlocationy()-y_average)*(particle[i].getdirectiony()-dy_average);
    z_sigma=z_sigma+(particle[i].getlocationz()-z_average)*(particle[i].getlocationz()-z_average);
    dz_sigma=dz_sigma+(particle[i].getdirectionz()-dz_average)*(particle[i].getdirectionz()-dz_average);
    zdz_sigma=zdz_sigma+(particle[i].getlocationz()-z_average)*(particle[i].getdirectionz()-dz_average);
  }
  x_sigma=x_sigma/particleLIVEnumber;
  dx_sigma=dx_sigma/particleLIVEnumber;
  xdx_sigma=xdx_sigma/particleLIVEnumber;
  correlationx=xdx_sigma/(sqrt(x_sigma*dx_sigma));
  y_sigma=y_sigma/particleLIVEnumber;
  dy_sigma=dy_sigma/particleLIVEnumber;
  ydy_sigma=ydy_sigma/particleLIVEnumber;
  correlationy=ydy_sigma/(sqrt(y_sigma*dy_sigma));


  z_sigma=z_sigma/particleLIVEnumber;
  dz_sigma=dz_sigma/particleLIVEnumber;
  zdz_sigma=zdz_sigma/particleLIVEnumber;
  emittancex=sqrt(x_sigma*dx_sigma-xdx_sigma*xdx_sigma);
  emittancey=sqrt(y_sigma*dy_sigma-ydy_sigma*ydy_sigma);
  emittancez=sqrt(z_sigma*dz_sigma-zdz_sigma*zdz_sigma);

  TBetax=x_sigma/emittancex;
  TBetay=y_sigma/emittancey;
  TBetaz=z_sigma/emittancez;
  TGammax=dx_sigma/emittancex;
  TGammay=dy_sigma/emittancey;
  TGammaz=dz_sigma/emittancez;
  if(TBetax*TGammax>1)
    TAlphax=sqrt(TBetax*TGammax-1);
  else
    TAlphax=0;
  if(TBetay*TGammay>1)
    TAlphay=sqrt(TBetay*TGammay-1);
  else
    TAlphay=0;
  if(TBetaz*TGammaz>1)
    TAlphaz=sqrt(TBetaz*TGammaz-1);
  else
    TAlphaz=0;

}
void BeamMap::output(char *p)
{
  ofstream out(p);
  for(int i=0;i<particleLIVEnumber;i++)
  {
    out<<particle[i].getlocationx()<<"	 ";
    out<<particle[i].getdirectionx()<<"	 ";
    out<<particle[i].getlocationy()<<"	 ";
    out<<particle[i].getdirectiony()<<"	 ";
    out<<particle[i].getlocationz()<<"	 ";
    out<<particle[i].getdirectionz()<<endl;
  }
  out.close();
}

//=========================================
double BeamMap::getsigmax()
{
  return x_sigma;
}
double BeamMap::getsigmay()
{
  return y_sigma;
}
double BeamMap::getsigmaz()
{
  return z_sigma;
}
double BeamMap::getsigmadx()
{
  return dx_sigma;
}
double BeamMap::getsigmady()
{
  return dy_sigma;
}
double BeamMap::getsigmadz()
{
  return dz_sigma;
}
double BeamMap::getsigmaxdx()
{
  return xdx_sigma;
}
double BeamMap::getsigmaydy()
{
  return ydy_sigma;
}
double BeamMap::getsigmazdz()
{
  return zdz_sigma;
}


double BeamMap::getTBetax()
{
  return TBetax;
}
double BeamMap::getTAlphax()
{
  return TAlphax;
}
double BeamMap::getTGammax()
{
  return TGammax;
}
double BeamMap::getTBetay()
{
  return TBetay;
}
double BeamMap::getTAlphay()
{
  return TAlphay;
}
double BeamMap::getTGammay()
{
  return TGammay;
}

double BeamMap::getTBetaz()
{
  return TBetaz;
}
double BeamMap::getTAlphaz()
{
  return TAlphaz;
}
double BeamMap::getTGammaz()
{
  return TGammaz;
}

double BeamMap::getbeta_beam_energy()
{
  return beta_beam_energy;
}
double BeamMap::getgamma_beam_energy()
{
  return gamma_beam_energy;
}
double BeamMap::getemittancex()
{
  return emittancex;
}
double BeamMap::getemittancey()
{
  return emittancey;
}
double BeamMap::getemittancez()
{
  return emittancez;
}
double BeamMap::getcorrelationx()
{
  return correlationx;
}
double BeamMap::getcorrelationy()
{
  return correlationy;
}
double BeamMap::getparticlenumber()
{
  return particleLIVEnumber;
}
void BeamMap::setparticlenumber(int a)
{
  particleLIVEnumber=a;
}

ParticleMap BeamMap::getparticle(int i)
{
  return particle[i];
}

void BeamMap::setparticle(int i,double a,double b,double c,double d,double e,double f)
{
  particle[i].setparticle(a,b,c,d,e,f);
}
void BeamMap::setlost(int i,int COUNT)
{
  swap(particle[i],particle[particleLIVEnumber-1]);
  particleLIVEnumber--;
  particle[particleLIVEnumber-1].lost=COUNT;
}



