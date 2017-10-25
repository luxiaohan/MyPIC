/*refer to trace3D RF Gap
 *Liuzc
 *Last Change @ 2015.3.25*/

#include "GapMap.h"
GapMap::~GapMap(){}
GapMap::GapMap(){}
GapMap::GapMap(double length,double E0T,double phase,double freq,double gamma)//(mm,MV/m,deg,Hz)
{
  this->length=length;
  this->phase=phase*M_PI/180;
  this->frequency=freq;
  this->wavelength=C_light/(freq)*1000;//mm
  this->E0T=E0T;
  this->Energygain=length/1000*E0T*1e6;
  if(phaseSYNflag==1)
    SliceMatrix=1;
  else
    SliceMatrix=100;
  setgamma(gamma);
  E0TFlag=1;
}

GapMap::GapMap(double length,const string &p,double phase_at_entry,double freq,double gamma)
{
  this->length=length;
  this->phase=phase_at_entry*M_PI/180;
  this->frequency=freq;
  this->wavelength=C_light/(freq)*1000;//mm
  F1.read(p);
  Slice=getlinenumber(p);
  if(phaseSYNflag==1)
    SliceMatrix=1;
  else
    SliceMatrix=100;
  setgamma(gamma);
  FieldFlag=1;
}

GapMap::GapMap(double tubelength,double gaplength,double radius,double phase_at_entry,double frequency,double gamma,int nharm)
{
/*
  this->length=length;
  this->phase=phase_at_entry*M_PI/180;
  this->frequency=frequency*1e6;
  this->wavelength=C_light/(frequency*1e6)*1000;//mm
  F1.read(p);
  Slice=getlinenumber(p);
  if(phaseSYNflag==1)
    SliceMatrix=1;
  else
    SliceMatrix=100;
  setgamma(gamma);
  FieldFlag=1;*/
}

void GapMap::Map()
{
  R=MatrixXd::Identity(6,6);
  DriftMap O1(length/SliceMatrix,gamma);
  for(int SliceCount=0;SliceCount<SliceMatrix;SliceCount++)
  {
    int h=1;						//harmonic number
    double bgi=sqrt(gamma*gamma-1);				//beta*gamma at begin
    double dw=Energygain/SliceMatrix*cos(phase);			//Enengy gain @MeV
    double gammaf=gamma+dw/(_MASS*C_light*C_light/abs(_Q));		//gamma at the end
    double bgf=sqrt(gammaf*gammaf-1);			//beta*gamma at end
    double gammaave=(gamma+gammaf)/2;
    double betaave=sqrt(1-1/gammaave/gammaave);
    //The changes in the normalized momentum components caused by the gap impulse
    double kxy=-1*M_PI*h*Energygain/SliceMatrix*sin(phase)	/	(_MASS*C_light*C_light/abs(_Q)*betaave*betaave*gammaave*gammaave*wavelength);

    /*
       double kdT_T=dphase/	(Energygain/SliceMatrix*_Q*sin(phase)/(_MASS*C_light*C_light*betaave*gammaave*gammaave));
       double kx=-abs(_Q)*Energygain/SliceMatrix*cos(phase)	/	(2*_MASS*C_light*C_light*betaave*betaave*gammaave*gammaave*gammaave)	*	(gammaave*gammaave+k_dT_T);
       double ky=-abs(_Q)*Energygain/SliceMatrix*cos(phase)	/	(2*_MASS*C_light*C_light*betaave*betaave*gammaave*gammaave*gammaave)	*	(gammazve*gammaave-k_dT_T);
     */
    double kz=2*M_PI*h*Energygain/SliceMatrix*sin(phase)/(_MASS*C_light*C_light/abs(_Q)*betaave*betaave*wavelength);
    MatrixXd RT=MatrixXd::Identity(6,6);
    /*
       cout<<bgi/bgf<<endl;
       cout<<bgi<<" aaa "<<bgf<<endl;
       cout<<gamma<<"  "<<gammaf<<endl;
       cout<<dw<<endl;
    */
    //xx
    RT(0,0)=1;
    RT(0,1)=0;
    RT(1,0)=kxy/bgf*1000;
    RT(1,1)=bgi/bgf;
    //yy
    RT(2,2)=1;
    RT(2,3)=0;
    RT(3,2)=kxy/bgf*1000;
    RT(3,3)=bgi/bgf;
    //zz
    RT(4,4)=1;
    RT(4,5)=0;
    RT(5,4)=kz/bgf*1000;
    RT(5,5)=bgi/bgf;
    R=RT*R;
    R=O1.getmap()*R;
    //cout<<phase<<"  "<<SliceCount<<"   "<<endl;
    phase+=2*M_PI*frequency/C_light*(length/SliceMatrix/1000)/beta;
    setgamma(gammaf);
  }
}
MatrixXd GapMap::getmap()
{
  return R;
}
MatrixXd GapMap::getmap(double x,double dx,double y,double dy,double z,double dp)
{
  double dphase=-2*M_PI/beta/wavelength*z;	//*sqrt(1+dx*dx/1e6/4+dy*dy/1e6/4)
  if(!(dphase<1e9))
  {
    cout<<"dphase error  "<<wavelength<<"  "<<beta<<"  "<<dphase<<endl;
    return R;
  }
  if(0)
    return R;
  //cout<<beta<<"  "<<dphase<<"  "<<wavelength<<"  "<<z<<endl;
  double phasetemp=phase;
  phase+=dphase;
  MatrixXd RTL=R;
  Map();
  R1=R;
  R=RTL;
  phase=phasetemp;
  return R1;
}

double GapMap::getEf()
{
  double dE,tempEf=1.0;
  vector<double>  Ef;
  vector<double> phaseT;
  double betaavg,gammaavg;
  double gammainitial=gamma,phaseinitial=phase;
  double gammatemp;
  double length_onestep=length/Slice/1;    	//mm to mm
  int be;
  for(int ll=0;ll<100;ll++)
  {
    Ef.push_back(tempEf);
    phaseT.push_back(0);
    setgamma(gammainitial);
    phase=phaseinitial;
    MatrixXd RL=MatrixXd::Identity(6,6);
    double phaseMatrix=0;
    for(int step=0;step<Slice;step++)
    {
      GapMap G1(length_onestep,F1.getfield(step)*Ef.back(),phase,frequency,gamma);
      G1.Map();
      RL=RL*G1.getmap();

      dE =_Q*F1.getfield(step)*1e6*Ef.back() * (length_onestep)  *  cos(phase);
      gammatemp=gamma+dE/(_MASS*C_light*C_light)/2;
      if(abs(gammatemp/gamma)>1.1)
      {
        be=abs(gammatemp*2/gamma);
       // be=1;
        length_onestep=length_onestep/be;
      }
      else
      {
        be=1;
      }
      //cout<<"  BE  "<<be<<endl;
      for(int lstep=0;lstep<be;lstep++)
      {
        dE =_Q*F1.getfield(step)*1e6*Ef.back() * (length_onestep)  *  cos(phase);
        gammatemp=gamma+dE/(_MASS*C_light*C_light)/2;
        setgamma(gammatemp);
        betaavg=beta;
        gammaavg=gamma;
        gamma+=dE/(_MASS*C_light*C_light)/2;
        setgamma(gamma);
        k=sqrt(abs(2*M_PI*_Q*F1.getfield(step)*1e6*Ef.back() * cos(phase) *  sin(phase)  /  (_MASS*C_light*C_light*pow(betaavg,3)*pow(gammaavg,3)*(wavelength/1))));
        phase+=2*M_PI  *  (length_onestep/1000)  /  (betaavg*C_light)  *  frequency;
        phaseT.back()+=k*length_onestep;
        if(1)
        {
          //cout<<step<<"  "<<phaseT.back()<<"  "<<RL(4,4)+RL(5,5)<<"  "<<phaseMatrix<<"  "<<phase<<endl;
        }
      }
      length_onestep=length_onestep*be;
    }
    phaseMatrix=acos((RL(4,4)+RL(5,5))/2);
    double alphaMT=(RL(4,4)-RL(5,5))/2/sin(phaseMatrix);
    double betaMT=RL(4,5)/sin(phaseMatrix);
    double gammaMT= (1+alphaMT*alphaMT)/betaMT;
    
    Vector3d tw;
    double PhaseAdvz=0,PhaseNowz=0;
    tw<<betaMT,alphaMT,gammaMT;
    for(int step=0;step<Slice;step++)
    {
      GapMap G1(length_onestep,F1.getfield(step)*Ef.back(),phase,frequency,gamma);
      G1.Map();
      tw=twissmapz(G1.getmap())*tw;
      PhaseAdvz=1/tw(0)*length_onestep/1000;	//tw(0)'s unit is mm/1000,not mm from Matrix
      //cout<<twissmapz(G1.getmap())<<endl;
      PhaseNowz+=PhaseAdvz;
    }
    
    
    cout<<Ef.back()<<" EF  "<<phaseT.back()<<"   "<<phaseMatrix<<"  "<<gamma<<"  "<<phase<<"  "<<PhaseNowz<<endl;
    //phaseT.back()=phaseMatrix;
    tempEf=Ef.back()/phaseT.back()*M_PI/180*89;
    if(!(phaseT.back()>=M_PI/180*90||phaseT.back()<M_PI/180*88))
      break;
    if(!(Ef.back()<1e9))
      break;
  }
  double temp;
  for(int i=0;i<Ef.size();i++) 
    for (int j =0;j< Ef.size() - 1 - i; j++) 
      if (Ef[j] < Ef[j + 1])
      {
        temp=Ef[j];
        Ef[j]=Ef[j+1];
        Ef[j+1]=temp;
        temp=phaseT[j];
        phaseT[j]=phaseT[j+1];
        phaseT[j+1]=temp;
      }
  for(int i=0;i<Ef.size();i++) 
    if(phaseT[i]<M_PI/180*90)
    {
      cout<<Ef[i]<<"  "<<phaseT[i]<<"  "<<i<<endl;
      return Ef[i];
    }
  return Ef[0];
}

/*
   void GapMap::FieldPromote();
   {
   }
 */
void GapMap::FieldMethod(double &x,double &px,double &y,double &py,double &z,double &pz)
{
  VectorXd XYZ(6);
  double pxt,pyt,pzt;
  double Ex,Ey,Ez,dw,dws,gamma0,beta0;
  double p1,p2,w1,w2;
  double dphase=-2*M_PI/beta/wavelength*z;

  double phasetemp=phase;
  double Energygaintemp=Energygain;
  MatrixXd RTL=R;
  length=length/Slice;

  XYZ<<x,px,y,py,z,pz;
  phase+=dphase;
  for(int FieldCount=0;FieldCount<Slice;FieldCount++)
  {
    E0T=F1.getfield(FieldCount);
    Energygain=length/1000*E0T*1e6;
    //cout<<E0T<<endl;   //to memory 1e1
    Map();
    XYZ=R*XYZ;
    dws=abs(_Q)*Energygain*length*cos(phase);
    gamma0=gamma;
    beta0=beta;
    gamma+=dws/(_MASS*C_light*C_light);
    setgamma(gamma);
    phase+=2*M_PI  *  length  /  ((beta0+beta)/2*C_light)  *  frequency;
  }

  x=XYZ(0);
  px=XYZ(1);
  y=XYZ(2);
  py=XYZ(3);
  z=XYZ(4);
  pz=XYZ(5);

  R=RTL;
  phase=phasetemp;
  Energygain=Energygaintemp;
  length=length*Slice;


  /*
     for(int step=0;step<Slice;step++)
     {
     Ex=0;
     Ey=0;
     Ez=F1.getfield(step)*cos(phase+dphase);
     dw=abs(_Q)*F1.getfield(step)*cos(phase+dphase)*length/Slice;
     dws=abs(_Q)*F1.getfield(step)*cos(phase)*length/Slice;
     gamma0=gamma;
     beta0=beta;
     gamma+=dws/(_MASS*C_light*C_light);
     setgamma(gamma);
     pxt=(px*beta0*C_light/1000+_Q*Ex/_MASS*length/Slice/(beta0*C_light))/(beta*C_light)*1000;
  //cout<<px<<"  PX  "<<pxt<<endl;
  pyt=(py*beta0*C_light/1000+_Q*Ey/_MASS*length/Slice/(beta0*C_light))/(beta*C_light)*1000;
  p1=(pz/1000+1)*gamma0*beta0*_MASS*C_light;
  w1=sqrt(p1*C_light*p1*C_light+(_MASS*C_light*C_light)*(_MASS*C_light*C_light));
  w2=w1+dw;
  p2=sqrt(w2*w2-(_MASS*C_light*C_light)*(_MASS*C_light*C_light))/C_light;
  pzt=(p2/(gamma*beta*_MASS*C_light)-1)*1000;
  x=x+(px+pxt)/1000/2*length/Slice;
  y=y+(py+pyt)/1000/2*length/Slice;
  z=z+(pz+pzt)/1000/2*length/Slice/gamma/gamma;
  px=pxt;
  py=pyt;
  pz=pzt;
  phase+=2*M_PI*length/Slice/wavelength;
  }
   */
}

void GapMap::print(){
  cout<<"The transfer matrix of Gap(Energy gain(EoTL) = "<<Energygain<<" V, phase = "<<phase/M_PI*180<<" deg,wavelength = "<<wavelength<<" ) is\n"<<R<<endl;
}

/*Funneling Gap===========================
  betaz=sqrt((1-1/gamma/gamma)/(1+dxp*dxp+dyp*dyp));
  k=q*Energygain*cos(phase);
//xx
R(0,0)=1;
R(0,1)=0;
R(1,0)=0;
R(1,1)=1;
//yy
R(2,2)=1;
R(2,3)=0;
R(3,2)=0;
R(3,3)=1;
//zz
R(4,4)=1;
R(4,5)=0;
R(5,4)=0;
R(5,5)=1;
//xz
R(0,4)=0;
R(0,5)=0;
R(1,4)=2*M_PI*abs(_Q)*Energygain*sin(phase)/(gamma*pow(beta,3)*wavelength*_MASS*C_light*C_light)*1000;
R(1,5)=0;
========================================*/

void GapMap::MakeMap(double tubelength,double gaplength,double radius,double frequency,int nharm)
{
  //transcode from LICHAO's fortran code
  int nn=1;
  double x,y,z,r;
  double lamada,tube,gap,beta,freq;
  double z0,z1,period,zgap,Egap;
  vector<double> coeff;
  vector<double> coeff1;
  vector<double> poten;
  vector<double> E,E1,Etemp;
  vector<double> W1,WR,AZ,AR,WZ;
  vector<double> EZ,ER,Btheta;
  vector<double> EEZ,EER,BBtheta;
  double Ex,Ey,Ezf,Erf,Bx,By,Bf,Bz;
  
  poten.resize(nn+1);
  tube=tubelength;
  gap=gaplength;
  freq=frequency;

  
  for(int i=0;i<poten.size();i++)
  {
    poten[i]=0;
  }
    poten[1]=1;
  E.resize(nharm+1);
  E1.resize(nharm*nn+1);
  coeff.resize(nharm*1);
  coeff1.resize(nn*nharm+1);
  Etemp.resize(nharm*nn+1);
  W1.resize(nharm*nn+1);
  WR.resize(nharm*nn+1);
  AZ.resize(nharm*nn+1);
  AR.resize(nharm*nn+1);
  WZ.resize(nharm*nn+1);
  EZ.resize(nharm*nn+1);
  ER.resize(nharm*nn+1);
  Btheta.resize(nharm*nn+1);
  EEZ.resize(nn+1);
  EER.resize(nn+1);
  BBtheta.resize(nn+1);
  
  z0=0.5*tube;
  z1=tube+gap+0.5*tube;
  zgap=tube+gap/2;
  period=tube+tube+gap*2;
  lamada=C_light/freq*1000;		//mm
  
  ofstream trackout("trackformat.dat");
  ofstream CSTout("CSTformat.dat");
  
  ULSCF(lamada,tube,gap,radius,poten,nharm,coeff,coeff1);
  for(int j=1;j<=nharm;j++)
  {
    E[j]=2*M_PI*j*coeff[j]/period;
    Egap+=E[j]*sin(2*M_PI*j*(zgap-z0)/period);
  }
  cout<<"Egap =  "<<Egap<<endl;
  
  ofstream coe("coeff.dat");
  for(int i=1;i<=nn;i++)
  {
    for(int j=1;j<=nharm;j++)
    {
      Etemp[(i-1)*nharm+j]=2*M_PI*j*coeff1[(i-1)*nharm+j]/period;
      coe<<setw(10)<<Etemp[(i-1)*nharm+j]<<"   ";
    }
    coe<<endl;
  }
  if(0)
  {
  int Slicex=25;
  int Slicey=25;
  int Slicez=201;
  for(int ix=1;ix<=Slicex;ix++)
  {
    x=-radius+2*radius/(Slicex-1)*(ix-1);
    for(int iy=1;iy<=Slicey;iy++)
    {
      y=-radius+2*radius/(Slicey-1)*(iy-1);
      r=sqrt(x*x+y*y);
      for(int iz=1;iz<=Slicez;iz++)
      {
        z=period/(Slicez-1)*(iz-1)/2;
        for(int i=0;i<EEZ.size();i++)
        {
          EEZ[i]=0;
          EER[i]=0;
          BBtheta[i]=0;
        }
        Ezf=0;
        Erf=0;
        Ex=0;
        Ey=0;
        Bx=0;
        By=0;
        Bz=0;
        Bf=0;
        for(int i=1;i<=nn;i++)
        {
          for(int j=1;j<=nharm;j++)
          {
            W1[(i-1)*nharm+j]=period/(j*lamada/i);
            E1[(i-1)*nharm+j]=Etemp[(i-1)*nharm+j];
            WR[(i-1)*nharm+j]=2.0*j*M_PI*r/period*sqrt(1-pow(W1[(i-1)*nharm+j],2));
            WZ[(i-1)*nharm+j]=2.0*j*M_PI/period*z*i;
            AZ[(i-1)*nharm+j]=bessi0(WR[(i-1)*nharm+j]);
	    AR[(i-1)*nharm+j]=bessi1(WR[(i-1)*nharm+j]);
	    EZ[(i-1)*nharm+j]=AZ[(i-1)*nharm+j]*Etemp[(i-1)*nharm+j]*sin(WZ[(i-1)*nharm+j]);
	    ER[(i-1)*nharm+j]=-AR[(i-1)*nharm+j]*Etemp[(i-1)*nharm+j]*cos(WZ[(i-1)*nharm+j])/sqrt(1-pow(W1[(i-1)*nharm+j],2));
	    Btheta[(i-1)*nharm+j]=-AR[(i-1)*nharm+j]*Etemp[(i-1)*nharm+j]*sin(WZ[(i-1)*nn+j])*W1[(i-1)*nharm+j]/sqrt(1-pow(W1[(i-1)*nharm+j],2))/3.0e+03;
	    EEZ[i]+=EZ[(i-1)*nharm+j];
	    EER[i]+=ER[(i-1)*nharm+j];
	    BBtheta[i]+=Btheta[(i-1)*nharm+j];
          }
          Ezf+=EEZ[i];
          Erf+=EER[i];
          Bf+=BBtheta[i];
        }
        if(r<1e-10)
        {
          Ex=0;
          Ey=0;
          Bx=0;
          By=0;
        }
        else
        {
          Ex=Erf*x/r;
          Ey=Erf*y/r;
          Bx=-Bf*y/r;
          By=Bf*x/r;
        }
        trackout<<setw(10)<<Ex<<"  "<<setw(10)<<Ey<<"  "<<setw(10)<<Ezf<<endl;
        CSTout<<setw(10)<<x<<"  "<<setw(10)<<y<<"  "<<setw(10)<<z<<"  "<<setw(10)<<Ex<<"  "<<setw(10)<<Ey<<"  "<<setw(10)<<Ezf<<"  "<<setw(3)<<0<<"  "<<setw(3)<<0<<"  "<<setw(3)<<0<<"  "<<endl;
      }
    }
  }
  }
  else
  {
    int Slicer=1;
    int Slicez=1001;
    for(int ix=1;ix<=Slicer;ix++)
    {
      r=sqrt(2.0)*0.9*radius/Slicer*(ix-1);
      for(int iz=1;iz<=Slicez;iz++)
      {
        z=period/(Slicez-1)*(iz-1)/2;
        for(int i=0;i<EEZ.size();i++)
        {
          EEZ[i]=0;
          EER[i]=0;
          BBtheta[i]=0;
        }
        Ezf=0;
        Erf=0;
        Ex=0;
        Ey=0;
        Bx=0;
        By=0;
        Bz=0;
        Bf=0;
        for(int i=1;i<=nn;i++)
        {
          for(int j=1;j<=nharm;j++)
          {
            W1[(i-1)*nharm+j]=period/(j*lamada/i);
            E1[(i-1)*nharm+j]=Etemp[(i-1)*nharm+j];
            WR[(i-1)*nharm+j]=2.0*j*M_PI*r/period*sqrt(1-pow(W1[(i-1)*nharm+j],2));
            WZ[(i-1)*nharm+j]=2.0*j*M_PI/period*z*i;
            AZ[(i-1)*nharm+j]=bessi0(WR[(i-1)*nharm+j]);
	    AR[(i-1)*nharm+j]=bessi1(WR[(i-1)*nharm+j]);
	    EZ[(i-1)*nharm+j]=AZ[(i-1)*nharm+j]*Etemp[(i-1)*nharm+j]*sin(WZ[(i-1)*nharm+j]);
	    ER[(i-1)*nharm+j]=-AR[(i-1)*nharm+j]*Etemp[(i-1)*nharm+j]*cos(WZ[(i-1)*nharm+j])/sqrt(1-pow(W1[(i-1)*nharm+j],2));
	    Btheta[(i-1)*nharm+j]=-AR[(i-1)*nharm+j]*Etemp[(i-1)*nharm+j]*sin(WZ[(i-1)*nn+j])*W1[(i-1)*nharm+j]/sqrt(1-pow(W1[(i-1)*nharm+j],2))/3.0e+03;
	    EEZ[i]+=EZ[(i-1)*nharm+j];
	    EER[i]+=ER[(i-1)*nharm+j];
	    BBtheta[i]+=Btheta[(i-1)*nharm+j];
          }
          Ezf+=EEZ[i];
          Erf+=EER[i];
          Bf+=BBtheta[i];
        }
        bool ONLY_VALUE=1;
        if(ONLY_VALUE=0)
          trackout<<setw(10)<<r<<"  "<<setw(10)<<z<<"  "<<setw(10)<<Erf<<"  "<<setw(10)<<Ezf<<"  "<<setw(10)<<Bf<<"  "<<endl;
        else
          trackout<<setw(10)<<Ezf<<"  "<<endl;
        if(0)
        {
          F1.addField(Ezf);
          Slice=Slicez;
        }
      }
    }
  }
}

void GapMap::ULSCF(double lamada,double tube,double gap,double radius,vector<double> poten,int nharm,vector<double> &coeff,vector<double> &coeff1)
{
  int nn=1;
  double temp,period,arg1,arg2,sq;
  period=tube*2+gap*2;
  for(int i=1;i<coeff.size();i++)	coeff[i]=0;
  for(int i=1;i<coeff1.size();i++)	coeff1[i]=0;
  for(int i=1;i<=nharm;i++)
  {
    arg1=i*M_PI*(tube+gap)/period;
    arg2=i*M_PI*(gap)/period;
    for(int j=1;j<=nn;j++)
    {
      temp=i;
      sq=sqrt(1-pow(period/(i*lamada/j),2));
      coeff[i]=coeff[i]+2*poten[j]/bessi0(2*M_PI*temp*radius*sq/period)*(tube+gap)/period*sin(arg1)/arg1*sin(arg2)/arg2;
      coeff1[(j-1)*nharm+i]=2*poten[j]/bessi0(2*M_PI*temp*radius*sq/period)*(tube+gap)/period*sin(arg1)/arg1*sin(arg2)/arg2;
      //cout<<lamada<<"  "<<coeff[i]<<"   "<<coeff1[(j-1)*nharm+i]<<"   "<<endl;
    }
  }
}
