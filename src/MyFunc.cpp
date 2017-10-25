#include "MyFunc.h"
#include <fstream>
#include <iostream>
#include <iomanip>

#include <iterator>
#include <cstdio> 
#include <cmath>
using namespace std;
using namespace Eigen;


const double min(const double *num,const int n)
{
  double x=num[0];
  for(int i=0;i<n;i++)
  {
    if(x>num[i]) x=num[i];
  }
  return x;
}

const double max(const double *num,const int n)
{
  double x=num[0];
  for(int i=0;i<n;i++)
  {
    if(x<num[i]) x=num[i];
  }
  return x;
}

const double ave(const double *num,const int n)
{
  double sum=0;
  for(int i=0;i<n;i++)
    sum+=num[i];
  return sum/n;
}

void StringSplit(const string &s,vector<string> &vec)
{
  if(vec.size()>0)
    vec.clear();
  int length = s.length();
  int start=0;
  char splitchar=' ';
  char splitchar2 = '	';
  for(int i=0;i<length;i++)
  {
    if (s[i]=='!')
    {i=length;}
    else if((s[i] == splitchar || s[i] == splitchar2) && i==start)
    {
      start++;
      //cout<<1<<endl;
    }
    else if(s[i] == splitchar || s[i] == splitchar2)
    {
      vec.push_back(s.substr(start,i - start));
      start = i+1;
      //cout<<2<<endl;
    }
    else if(i == length-1)
    {
      vec.push_back(s.substr(start,i+1 - start));
      //    cout<<3<<endl;
    }
  }
}

int getlinenumber(const string &p)
{
  //	fstream fsCopee(p, ios::binary | ios::in ) ;
  //	fstream fsCoper("temp.dat", ios::binary | ios::out ) ;
  //	fsCoper << fsCopee.rdbuf() ;
  ifstream file(p);
  string str;
  int count = 0;
  while (file.peek()!=EOF) {
    getline(file, str);
    remove(str.begin(), str.end(), ' ');
    remove(str.begin(), str.end(), '\t');
    if (str.length() > 0) 
    {
      count ++;
    }
  }
  //cout<<"particle number is "<<count<<endl;
  file.close();
  return count;
}

double gaussrand()//
{
  static double V1, V2, S;
  static int phase = 0;
  double X;

  if ( phase == 0 ) {
    do {
      double U1 = (double)rand() / RAND_MAX;
      double U2 = (double)rand() / RAND_MAX;

      V1 = 2 * U1 - 1;
      V2 = 2 * U2 - 1;
      S = V1 * V1 + V2 * V2;

    } while(S >= 1 || S == 0);

    X = V1 * sqrt(-2 * log(S) / S);
  } else
    X = V2 * sqrt(-2 * log(S) / S);

  phase = 1 - phase;

  return X;
}

double bessi0(double x)
{

  double result=0;
  double nx=1;
  for(int i=0;i<20;++i)
  {
    if(i>0)
    {
      nx/=i;
    }
    result+=pow(x/2,2*i)/(tgamma(i+1))*nx;
  }
  return result;
  /*
  double X1,bess;
  if((x=abs(x))<=3.75)
  {
    X1=pow(x/3.75,2);
    bess=1.+X1*(3.5156229+X1*(3.0899424+X1*(1.2067492+X1*(0.2659732+X1*(0.0360768+X1*0.0045813)))));
  }
  else
  {
    X1=3.75/x;
    bess=(0.398942280+X1*(0.013285917+X1*(0.002253187-X1*(0.001575649-X1*(0.009162808-X1*(0.020577063-X1*(0.026355372-X1*(0.016476329-X1*0.003923767))))))))*exp(x)/sqrt(x);
  }
  return bess;
  */
}

double bessi1(double x)
{

  double result=0;
  double nx=1;
  for(int i=0;i<20;++i)
  {
    if(i>0)
    {
      nx/=i;
    }
    result+=pow(x/2,2*i+1)/(tgamma(i+1+1))*nx;
  }
  return result;
  /*
  double ax,ans,y;
  if((ax=fabs(x)) <= 3.75)
  {
    y=x/3.75;
    y*=y;
    ans=ax*(0.5+y*(0.87890594+y*(0.51498869+y*(0.15084934
              +y*(0.2658733e-1+y*(0.301532e-2+y*0.32411e-3))))));
  }
  else
  {
    y=3.75/ax;
    ans=0.2282967e-1+y*(-0.2895312e-1+y*(0.1787654e-1-y*0.420059e-2));
    ans=0.39894228+y*(-0.3988024e-1+y*(-0.362018e-2+y*(0.163801e-2+y*(-0.1031555e-1+y*ans))));
    ans *= (exp(ax)/sqrt(ax));
  }
  return x < 0.0 ? -ans : ans;
  */
}

double Bessin(int n,double x)
{
  double result=0;
  double nx=1;
  for(int i=0;i<20;++i)
  {
    if(i>0)
    {
      nx/=i;
    }
    result+=pow(x/2,2*i+n)/(tgamma(i+1+n))*nx;
  }
  return result;
}

MatrixXd twissmapx(const MatrixXd &R)//beta alpha gamma,in x direction test
{
  MatrixXd RTEMP=MatrixXd::Zero(3,3);
  RTEMP<<R(0,0)*R(0,0),-2*R(0,0)*R(0,1),R(0,1)*R(0,1),
    -R(0,0)*R(1,0),R(0,0)*R(1,1)+R(0,1)*R(1,0),-R(0,1)*R(1,1),
    R(1,0)*R(1,0),-2*R(1,0)*R(1,1),R(1,1)*R(1,1);
  return RTEMP;
}
MatrixXd twissmapy(const MatrixXd &R)//beta alpha gamma,in y direction test
{
  MatrixXd RTEMP=MatrixXd::Zero(3,3);
  RTEMP<<R(2,2)*R(2,2),-2*R(2,2)*R(2,3),R(2,3)*R(2,3),
    -R(2,2)*R(3,2),R(2,2)*R(3,3)+R(2,3)*R(3,2),-R(2,3)*R(3,3),
    R(3,2)*R(3,2),-2*R(3,2)*R(3,3),R(3,3)*R(3,3);
  return RTEMP;
}
MatrixXd twissmapz(const MatrixXd &R)//beta alpha gamma,in z direction test
{
  MatrixXd RTEMP=MatrixXd::Zero(3,3);
  RTEMP<<R(4,4)*R(4,4),-2*R(4,4)*R(4,5),R(4,5)*R(4,5),
    -R(4,4)*R(5,4),R(4,4)*R(5,5)+R(4,5)*R(5,4),-R(4,5)*R(5,5),
    R(5,4)*R(5,4),-2*R(5,4)*R(5,5),R(5,5)*R(5,5);
  return RTEMP;
}

void printvector(ostream &par,const vector<vector<double>  > &pv)
{
  //ostream_iterator<double> outpar(par,"\t");
  //for_each(pv.begin(),pv.end(), [&par,&outpar](const vector<double> &a){  copy(a.begin(),a.end(),outpar);par<<endl;});
  //for(const auto &a:pv) {copy(a.begin(),a.end(),outpar);  par<<endl;}
  par.setf(ios::scientific, ios::floatfield);
  for(int i =0;i<pv.size();++i)
  {
    for(int j=0;j<pv[i].size();++j)
    {
      par << setprecision(8) <<setiosflags(ios::showpos)<<pv[i][j]<<"   ";
    }
    par <<endl;
  }
}

