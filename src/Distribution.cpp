#include "Distribution.h"
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cstdio> 
#include <cstdlib>
#include <cmath>
#include <complex>
#include <Eigen/Dense>
#include <gsl/gsl_sf_lambert.h>

#include "MyFunc.h"
using namespace std;

void CKV_Beam(const vector<double> &param, int particlenumber)
{}
void KV2D_T(const vector<double> &twiss,vector<vector<double> > &xy)
{  
  //vector twiss:
  //Emitx, Betax, Alphax, 
  //Emity, Betay, Alphay, 
  if(twiss.size()<6)    return;

  double sigx	= twiss.at(1)   *	twiss.at(0)			/4;
  double sigpx	= ( pow(twiss.at(2),2)+1 )/twiss.at(1)*twiss.at(0)	/4;
  double sigxpx	= -twiss.at(2)*twiss.at(0)				/4;
  double sigy	= twiss.at(4) 	*	twiss.at(3)			/4;
  double sigpy	= ( pow(twiss.at(5),2)+1 )/twiss.at(4)*twiss.at(3)	/4;
  double sigypy	= -twiss.at(5)*twiss.at(3)				/4;
  //  cout<<"sig : "<<sigx<<"  "<<sigpx<<"  "<<sigxpx<<endl;
  KV2D_S(sigx,sigpx,sigxpx,sigy,sigpy,sigypy,xy);
}


void KV3D_T(const vector<double> &twiss,vector<vector<double> > &xy)
{
  //vector twiss:
  //Emitx, Betax, Alphax, 
  //Emity, Betay, Alphay, 
  //Emitz, Betaz, Alphaz	in (z,z') plane
  if(twiss.size()<9)    
  {
    cout<<"KV3D_T input parameter error!"<<endl;
    return;
  }
  KV2D_T(twiss,xy);
  double sigz	= twiss.at(7)	*	twiss.at(6)			/4;
  double sigpz	= ( pow(twiss.at(8),2)+1 )/twiss.at(7)*twiss.at(6)	/4;
  double sigzpz	= twiss.at(8)*twiss.at(6)				/4;
  vector<double> twissz(6);
  vector<vector<double> > z=xy;
  for(int i=0;i<6;++i)
  {
    if(i<3)
    {
      twissz[i]=twiss[i+6];
      //cout<<twissz[i]<<endl;
    }
    else
    {
      twissz[i]=twiss[i+3];
    }
  }
  KV2D_T(twissz,z);

  for(int i=0;i<xy.size();++i)
  {
    xy[i][4]=z[i][0];
    xy[i][5]=z[i][1];
  }

}


void KV2D_S(double sig_x,double sig_px,double sig_xpx,double sig_y,double sig_py,double sig_ypy,vector<vector<double> > &xy)
{
  if(sig_x*sig_px<=sig_xpx*sig_xpx||sig_y*sig_py<=sig_ypy*sig_ypy)
  {
    cout<<"ERROR!!!   @sigma!"<<endl;
    return;
  }
  Eigen::MatrixXd CM = Eigen::MatrixXd::Zero(2,2);
  Eigen::VectorXd eivals;

  double x1=2*sig_xpx/sqrt(sig_px);
  double px1=2*sqrt(sig_px);
  double x2=2*sqrt(sig_x);
  double px2=2*sig_xpx/sqrt(sig_x);
  double Bx=-2*px2/px1/(x2*px1-x1*px2);
  double Cx=x2/px1/(x2*px1-x1*px2);
  double Ax=(1+px2*px2*Cx)/x2/x2;
  CM<<Ax,   Bx/2,
    Bx/2, Cx;
  eivals = CM.selfadjointView<Eigen::Lower>().eigenvalues();
  if(sig_x<sig_px)
  {
    double temp = eivals(0);
    eivals(0)=eivals(1);
    eivals(1)=temp;
  }
  double ax=sqrt(abs(1/eivals(0)));
  double bx=sqrt(abs(1/eivals(1)));
  double phix;
  if(Ax!=Cx)
    phix=1.0/2*atan(Bx/(Ax-Cx))*180/M_PI;
  else
    phix=0;
  double emix=ax*bx;

  double y1=2*sig_ypy/sqrt(sig_py);
  double py1=2*sqrt(sig_py);
  double y2=2*sqrt(sig_y);
  double py2=2*sig_ypy/sqrt(sig_y);
  double By=-2*py2/py1/(y2*py1-y1*py2);
  double Cy=y2/py1/(y2*py1-y1*py2);
  double Ay=(1+py2*py2*Cy)/y2/y2;
  CM<<Ay,   By/2,
    By/2, Cy;
  eivals = CM.selfadjointView<Eigen::Lower>().eigenvalues();
  if(sig_y<sig_py)
  {
    double temp = eivals(0);
    eivals(0)=eivals(1);
    eivals(1)=temp;
  }
  double ay=sqrt(1/eivals(0));
  double by=sqrt(1/eivals(1));
  double phiy;
  if(Ay!=Cy)
    phiy=1.0/2*atan(By/(Ay-Cy))*180/M_PI;
  else
    phiy=0;
  double emiy=ay*by;
  /*
     cout<<"x1 = "<<x1<<" , px1= "<<px1<<endl;
     cout<<"x2 = "<<x2<<" , px2= "<<px2<<endl;
     cout<<Ax<<"  "<<Bx<<"  "<<Cx<<"  ABC"<<endl;
     cout<<phix<<"  phi  "<<phiy<<endl;
     cout<<emix<<"  "<<ax<<"  "<<emiy<<"  "<<ay<<"  "<<phix<<"  "<<phiy<<"  "<<"  "<<p<<endl;;
     cout<<"emit x exp : "<<4*sqrt(sig_x*sig_px-sig_xpx*sig_xpx)<<endl;
   */
  KV2D_A(vector<double>{emix,ax,bx,phix,emiy,ay,by,phiy},xy);
}

void KV2D_A(const vector<double> &param,vector<vector<double> > &xy)
{
  vector<double> paramRet=distribution_base1(param);
  double Factor=param.at(0);
  double G=0;
  for(auto num=xy.begin();num!=xy.end();++num)
  {
    //G=(double)rand() / RAND_MAX;
    paramRet.back()=Factor;
    distribution_base2(paramRet,*num);
  }
}

// ==============================
void KVdistributionT(const vector<double> &twiss,int particlenumber,char *p)
{
  //call
  vector<vector<double>  > particle(particlenumber, vector<double>(6));
  KV2D_T(twiss,particle);
  //print
  ofstream par(p);
  printvector(par,particle);
  par.close();
}


void KVdistributionS(double sig_x,double sig_px,double sig_xpx,double sig_y,double sig_py,double sig_ypy,int particlenumber,char *p)
{
  //call
  vector<vector<double>  > particle(particlenumber, vector<double>(6));
  KV2D_S(sig_x,sig_px,sig_xpx,sig_y,sig_py,sig_ypy,particle);
  //print
  ofstream par(p);
  printvector(par,particle);
  par.close();
}

void KVdistributionA(const vector<double> &param,int particlenumber,char *p)
{
// vector param consist the following sequencely:
//   emittence_X
//   semiaxis_X1
//   semiaxis_X2
//   tilt_X
//   emittence_Y
//   semiaxis_Y1
//   semiaxis_Y2
//   tilt_Y
  //call
  vector<vector<double>  > particle(particlenumber, vector<double>(6));
  KV2D_A(param,particle);
  //print
  ofstream par(p);
  printvector(par,particle);
  par.close();
}


void KVdistribution6T(const vector<double> &twiss,int particlenumber,char *p)
{
  //vector twiss:
  //Emitx, Betax, Alphax, 
  //Emity, Betay, Alphay, 
  //Emitz, Betaz, Alphaz
  if(twiss.size()<9)
    return;
  double sigx	= twiss.at(1)   *	twiss.at(0)			/4;
  double sigpx	= ( pow(twiss.at(2),2)+1 )/twiss.at(1)*twiss.at(0)	/4;
  double sigxpx	= twiss.at(2)*twiss.at(0)				/4;
  double sigy	= twiss.at(4) 	*	twiss.at(3)			/4;
  double sigpy	= ( pow(twiss.at(5),2)+1 )/twiss.at(4)*twiss.at(3)	/4;
  double sigypy	= twiss.at(5)*twiss.at(3)				/4;
  double sigz	= twiss.at(7)	*	twiss.at(6)			/4;
  double sigpz	= ( pow(twiss.at(8),2)+1 )/twiss.at(7)*twiss.at(6)	/4;
  double sigzpz	= twiss.at(8)*twiss.at(6)				/4;
  cout<<"sig : "<<sqrt(sigx)<<"  "<<sqrt(sigy)<<"  "<<sigxpx<<endl;
  KVdistribution6S(sigx,sigpx,sigxpx,sigy,sigpy,sigypy,sigz,sigpz,sigzpz,particlenumber,p);
}

void KVdistribution6S(double sig_x,double sig_px,double sig_xpx,double sig_y,double sig_py,double sig_ypy,double sig_z,double sig_pz,double sig_zpz,int particlenumber,char *p)
{
  if(sig_x*sig_px<sig_xpx*sig_xpx||sig_y*sig_py<sig_ypy*sig_ypy||sig_z*sig_pz<sig_zpz*sig_zpz)
  {
    cout<<"error sigma minu"<<endl;
    return;
  }
  Eigen::MatrixXd CM = Eigen::MatrixXd::Zero(2,2);
  //  Eigen::MatrixXd DM = Eigen::MatrixXd::Zero(3,3);
  Eigen::VectorXd eivals;
  double x1=2*sig_xpx/sqrt(sig_px);
  double px1=2*sqrt(sig_px);
  double x2=2*sqrt(sig_x);
  double px2=2*sig_xpx/sqrt(sig_x);
  //  double Ax=px2/x1/(x2*px1-x1*px2);
  double Bx=-2*px2/px1/(x2*px1-x1*px2);
  double Cx=x2/px1/(x2*px1-x1*px2);
  double Ax=(1+px2*px2*Cx)/x2/x2;
  CM<<Ax,   Bx/2,
    Bx/2, Cx;
  /*
     DM<<Ax,   Bx/2, 0,
     Bx/2, Cx  , 0,
     0  ,  0  ,  1;
   */
  eivals = CM.selfadjointView<Eigen::Lower>().eigenvalues();
  if(sig_x<sig_px)
  {
    double temp = eivals(0);
    eivals(0)=eivals(1);
    eivals(1)=temp;
  }
  double ax=sqrt(abs(1/eivals(0)));
  double bx=sqrt(abs(1/eivals(1)));
  double phix;
  if(Ax!=Cx)
    phix=1.0/2*atan(Bx/(Ax-Cx))*180/M_PI;
  else
    phix=0;
  double emix=ax*bx;
  
  double y1=2*sig_ypy/sqrt(sig_py);
  double py1=2*sqrt(sig_py);
  double y2=2*sqrt(sig_y);
  double py2=2*sig_ypy/sqrt(sig_y);
  double By=-2*py2/py1/(y2*py1-y1*py2);
  double Cy=y2/py1/(y2*py1-y1*py2);
  double Ay=(1+py2*py2*Cy)/y2/y2;
  CM<<Ay,   By/2,
    By/2, Cy;
  eivals = CM.selfadjointView<Eigen::Lower>().eigenvalues();
  if(sig_y<sig_py)
  {
    double temp = eivals(0);
    eivals(0)=eivals(1);
    eivals(1)=temp;
  }
  double ay=sqrt(1/eivals(0));
  double by=sqrt(1/eivals(1));
  double phiy;
  if(Ay!=Cy)
    phiy=1.0/2*atan(By/(Ay-Cy))*180/M_PI;
  else
    phiy=0;
  double emiy=ay*by;
  
  double z1=2*sig_zpz/sqrt(sig_pz);
  double pz1=2*sqrt(sig_pz);
  double z2=2*sqrt(sig_z);
  double pz2=2*sig_zpz/sqrt(sig_z);
  double Bz=-2*pz2/pz1/(z2*pz1-z1*pz2);
  double Cz=z2/pz1/(z2*pz1-z1*pz2);
  double Az=(1+pz2*pz2*Cz)/z2/z2;
  CM<<Az,   Bz/2,
    Bz/2, Cz;
  eivals = CM.selfadjointView<Eigen::Lower>().eigenvalues();
  if(sig_z<sig_pz)
  {
    double temp = eivals(0);
    eivals(0)=eivals(1);
    eivals(1)=temp;
  }
  double az=sqrt(1/eivals(0));
  double bz=sqrt(1/eivals(1));
  double phiz;
  if(Az!=Cz)
    phiz=1.0/2*atan(Bz/(Az-Cz))*180/M_PI;
  else
    phiz=0;
  double emiz=az*bz;
  
  /*
     cout<<"x1 = "<<x1<<" , px1= "<<px1<<endl;
     cout<<"x2 = "<<x2<<" , px2= "<<px2<<endl;
     cout<<Ax<<"  "<<Bx<<"  "<<Cx<<"  ABC"<<endl;
     cout<<phix<<"  phi  "<<phiy<<endl;
     cout<<emix<<"  "<<ax<<"  "<<emiy<<"  "<<ay<<"  "<<phix<<"  "<<phiy<<"  "<<particlenumber<<"  "<<p<<endl;;
     cout<<"emit x exp : "<<4*sqrt(sig_x*sig_px-sig_xpx*sig_xpx)<<endl;
   */
  KVdistribution6A(emix,ax,emiy,ay,emiz,az,phix,phiy,phiz,particlenumber,p);
  /*
     particleread(p);
     cout<<beam.getsigmax()<<endl;
     cout<<beam.getsigmadx()<<endl;
     cout<<beam.getsigmaxdx()<<endl;
   */
}

void KVdistribution6A(double emitx,double betax,double emity,double betay,double emitz,double betaz,double alpx,double alpy,double alpz,int particlenumber,char *p)
{
  ofstream par(p);
  double co1,co2;
  double ax,axs,ay,ays,az,azs,emitratio,emitratio2,F;
  double x,xs,y,ys,z,zs;
  double anglebetax,anglebetay,anglebetaz,zetax,zetay,zetaz;
  double alphax=alpx/180*M_PI;
  double alphay=alpy/180*M_PI;
  double alphaz=alpz/180*M_PI;
  double x1=betax;
  double x2=emitx/betax;
  double y1=betay;
  double y2=emity/betay;
  double z1=betaz;
  double z2=emitz/betaz;
  ax=sqrt(x1/x2*pow(cos(alphax),2)+x2/x1*pow(sin(alphax),2));
  axs=(1/ax/2)*((x1/x2)-(x2/x1))*sin(2*alphax);
  ay=sqrt(y1/y2*pow(cos(alphay),2)+y2/y1*pow(sin(alphay),2));
  ays=(1/ay/2)*(y1/y2-y2/y1)*sin(2*alphay);
  az=sqrt(z1/z2*pow(cos(alphaz),2)+y2/y1*pow(sin(alphaz),2));
  azs=(1/az/2)*(z1/z2-z2/z1)*sin(2*alphaz);
  emitx=x1*x2;
  emity=y1*y2;
  emitz=z1*z2;
  emitratio=emitx/emity;
  emitratio2=emitx/emitz;
  F=emitx;
  //  cout<<"F! "<<F<<endl;
  for(int num=0;num<particlenumber;num++)
  {
    co1=(double)rand() / RAND_MAX;
    co2=(double)rand() / RAND_MAX;
/*    if(co1>co2)
    {
      double t=co1;
      co1=co2;
      co2=t;
    }*/
    zetax=sqrt(co1*F);
    zetay=sqrt((1-co1)*F/emitratio);
    zetaz=sqrt((co2)*F/emitratio2);
    anglebetax=(double)rand() / RAND_MAX*2*M_PI;
    anglebetay=(double)rand() / RAND_MAX*2*M_PI;
    anglebetaz=(double)rand() / RAND_MAX*2*M_PI;
    x=zetax*ax*cos(anglebetax);
    xs=zetax*(axs*cos(anglebetax)-sin(anglebetax)/ax);
    y=zetay*ay*cos(anglebetay);
    ys=zetay*(ays*cos(anglebetay)-sin(anglebetay)/ay);
    z=zetaz*az*cos(anglebetaz);
    zs=zetaz*(azs*cos(anglebetaz)-sin(anglebetaz)/az);
    par<<setw(10)<<x<<"   "<<setw(10)<<xs<<"   "<<setw(10)<<y<<"   "<<setw(10)<<ys<<"   "<<setw(10)<<z<<"   "<<setw(10)<<zs<<endl;
  }
  par.close();
}

void distribution(double emitx,double betax,double emity,double betay,double rx,double ry,double rxy,double rxtyt,double rxyt,double rxty,int particlenumber,char *p)
{
  //x,x'   y,y'
  double paramRet=sqrt(emitx*betax);
  double b=sqrt(emitx/betax);
  double c=sqrt(emity*betay);
  double d=sqrt(emity/betay);
  ofstream par(p);
  ofstream par1("particletests.dat");
  double U1,U2,U3,U4,X1,X2,X3,X4;
  Eigen::MatrixXd C=Eigen::MatrixXd::Identity(4,4);
  Eigen::MatrixXd m = Eigen::MatrixXd::Identity(4,4);
  m(0,0)=paramRet*paramRet;
  m(1,1)=b*b;
  m(0,1)=rx*paramRet*b;
  m(1,0)=rx*paramRet*b;
  m(2,2)=c*c;
  m(3,3)=d*d;
  m(2,3)=ry*c*d;
  m(3,2)=ry*c*d;
  m(0,2)=rxy*paramRet*c;
  m(2,0)=rxy*paramRet*c;
  m(0,3)=rxyt*paramRet*d;
  m(3,0)=rxyt*paramRet*d;
  m(1,2)=rxty*b*c;
  m(2,1)=rxty*b*c;
  m(1,3)=rxtyt*b*d;
  m(3,1)=rxtyt*b*d;		//or m(3,1)=rxtyt*rxtyt*b*d;?
  C = m.llt().matrixL();
  //cout<<C<<endl;
  for(int number_of_particles=0;number_of_particles<particlenumber;number_of_particles++)//
  {
    double x[4]={0,0,0,0};
    double mu[4]={0,0,0,0};//:expectation
    double u[4];
    for(int i=0;i<4;i++)
    {
      u[i]=gaussrand();
      par1<<setw(10)<<u[i]<<"   ";
    }
    par1<<setw(10)<<0<<"   "<<setw(10)<<0<<endl;
    for(int k=0;k<4;k++)
    {
      for(int i=0;i<=k;i++)
      {
        x[k]=x[k]+C(k,i)*u[i];
      }
      x[k]=x[k]+mu[k];
      par<<setw(10)<<x[k]<<"   ";
    }
    double z=(double)rand() / RAND_MAX*1;
    double zs=gaussrand()*1;
    par<<setw(10)<<z<<"   "<<setw(10)<<zs<<endl;
  }
  par.close();
  par1.close();
}

void WB2D_A(const vector<double> &param,vector<vector<double> > &xy)
{
  vector<double> paramRet=distribution_base1(param);
  double Factor=param.at(0)/2*3;
  double G=0;
  for(auto num=xy.begin();num!=xy.end();++num)
  {
    G=(double)rand() / RAND_MAX;
    paramRet.back()=Factor*sqrt(G);
    distribution_base2(paramRet,*num);
  }
}

void PB2D_A(const vector<double> &param,vector<vector<double> > &xy)
{
  vector<double> paramRet=distribution_base1(param);
  double Factor=param.at(0)*2;
  double G=0;
  complex<double> temp2   =1.0 + complex<double>(0.0,1.0) * sqrt(3);
  complex<double> temp1   =1.0 - complex<double>(0.0,1.0) * sqrt(3);
  complex<double> temp3(0,0);
  for(auto num=xy.begin();num!=xy.end();++num)
  {
    G=(double)rand() / RAND_MAX;
    temp3=pow(1 - 2*G + 2.0*sqrt(complex<double>(-G + pow(G,2),0.0)),1.0/3);
    paramRet.back()=Factor*(0.5 + (temp1/(4.0*temp3)).real()  +  (temp2*temp3/4.0).real());
    distribution_base2(paramRet,*num);
  }
}

void GS2D_A(const vector<double> &param,vector<vector<double> > &xy)
{
// vector param consist the following sequencely:
//   emittence_X
//   semiaxis_X1
//   semiaxis_X2
//   tilt_X
//   emittence_Y
//   semiaxis_Y1
//   semiaxis_Y2
//   tilt_Y
  vector<double> paramRet=distribution_base1(param);
  double Factor=param.at(0)/2;
  double G=0;
  for(auto num=xy.begin();num!=xy.end();++num)
  {
    G=(double)rand() / RAND_MAX;
    paramRet.back()=Factor*(-1-gsl_sf_lambert_Wm1((G-1)/M_E));
    //cout<<G<<"  "<<Factor<<"  "<<paramRet.back()/Factor<<endl;
    distribution_base2(paramRet,*num);
  }
}


vector<double> distribution_base1(const vector<double> &param)
{
// vector param consist the following sequencely:
//   emittence_X
//   semiaxis_X1
//   semiaxis_X2
//   tilt_X
//   emittence_Y
//   semiaxis_Y1
//   semiaxis_Y2
//   tilt_Y
  double ax=0,axs=0,ay=0,ays=0;
  ax=sqrt(  param.at(1)/param.at(2)*pow(cos(param.at(3)*M_PI/180),2)+param.at(2)/param.at(1)*pow(sin(param.at(3)*M_PI/180),2)  );
  axs=1.0/(ax*2)*((param.at(1)/param.at(2))-(param.at(2)/param.at(1)))*sin(2*param.at(3)*M_PI/180);
  ay=sqrt(  param.at(5)/param.at(6)*pow(cos(param.at(7)*M_PI/180),2)+param.at(6)/param.at(5)*pow(sin(param.at(7)*M_PI/180),2)  );
  ays=1.0/(ay*2)*((param.at(5)/param.at(6))-(param.at(6)/param.at(5)))*sin(2*param.at(7)*M_PI/180);
  double emitratio=param.at(0)/param.at(4);
  return vector<double>{ax,axs,ay,ays,emitratio,0.0};//the last place is reserved for F 
}

void distribution_base2(const vector<double> &par,vector<double> &xy)
{
    double zetax=0,zetay=0,anglebetax=0,anglebetay=0;
    zetax=sqrt((double)rand() / RAND_MAX*par.back());
    zetay=sqrt((par.back()-zetax*zetax)/par.at(4));
    anglebetax=(double)rand() / RAND_MAX*2*M_PI;
    anglebetay=(double)rand() / RAND_MAX*2*M_PI;
    xy.at(0)=zetax  *   par.at(0)*cos(anglebetax);
    xy.at(1)=zetax  *  (par.at(1)*cos(anglebetax)-sin(anglebetax)/par.at(0));
    xy.at(2)=zetay  *   par.at(2)*cos(anglebetay);
    xy.at(3)=zetay  *  (par.at(3)*cos(anglebetay)-sin(anglebetay)/par.at(2));
    xy.at(4)=0.0;
    xy.at(5)=0.0;
}


//according to semi-axis length and tilt, get the coordinate of the point of maximum of x max y
//fun()
//{
//  double limity_x,limity_y,limitx_x,limitx_y;
//  if(alpx!=0&&alpx!=90)
//  {
//    limitx_x=abs(sqrt((x2*x2+x1*x1/(tan(alphax)*tan(alphax)))/(1+1/(tan(alphax)*tan(alphax)))));
//    limitx_y=abs(((1/tan(alphax)*sqrt((x2*x2+x1*x1/(tan(alphax)*tan(alphax))))/x2/x2/(1/x1/x1+1/x2/x2/(tan(alphax)*tan(alphax))))-abs(limitx_x*cos(alphax)))/sin(alphax));
//    limity_y=abs(sqrt((x2*x2+x1*x1*tan(alphax)*tan(alphax))/(1+tan(alphax)*tan(alphax))));
//    limity_x=abs(((tan(alphax)*sqrt((x2*x2+x1*x1*(tan(alphax)*tan(alphax))))/x2/x2/(1/x1/x1+1/x2/x2*(tan(alphax)*tan(alphax))))-abs(limity_y*sin(alphax)))/cos(alphax));
//  }
//  else if(alpx==0)
//  {
//    limitx_x=x1;
//    limitx_y=0;
//    limity_y=x2;
//    limity_x=0;
//  }
//  else
//  {
//    limitx_x=x2;
//    limitx_y=0;
//    limity_y=x1;
//    limity_x=0;
//  }
//  emitx=sqrt(limitx_x*limity_y*limitx_x*limity_y-limitx_x*limitx_y*limitx_x*limitx_y);
//  /*
//     if(myid==0)
//     {
//     cout<<"x1re = "<<limity_x<<" , px1re= "<<limity_y<<endl;
//     cout<<"x2re = "<<limitx_x<<" , px2re= "<<limitx_y<<endl;
//     cout<<limitx_x*limitx_y<<"  "<<limity_x*limity_y<<endl;
//     }
//   */
//}
