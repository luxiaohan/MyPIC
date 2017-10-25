#include <cmath>
#include <iostream>
#include <vector>

using std::vector;
using std::cout;
using std::endl;



void fourier(vector<double> &in,int sign)
{
  const int n =in.size();
  const int nn=n/2;
  int j=1,m,mmax,istep;
  double tempr,tempi;
  double theta,wpr,wpi,wr,wi;
  for(int i=1;i<=n;i+=2)
  {
    if(j>i)
    {
      tempr  = in[j-1];
      tempi  = in[j  ];
      in[j-1]= in[i-1];
      in[j  ]= in[i  ];
      in[i-1]= tempr;
      in[i  ]= tempi;
    }
    m=nn;
    while(m>=2 && j>m)
    {
      j-=m;
      m/=2;
    }
    j+=m;
  }
  mmax=2;
  while(n>mmax)
  {
    istep= 2*mmax;
    theta= 2*M_PI/(sign*mmax);
    wpr  =-2*pow( sin(theta/2) , 2);
    wpi  = sin(theta);
    wr   = 1;
    wi   = 0;
    for(m = 1; m<=mmax ; m+=2)
    {
      for(int i=m; i<=n;i+=istep)
      {
	j=i+mmax;
	tempr=wr*in[j-1] - wi*in[j  ];
	tempi=wr*in[j  ] + wi*in[j-1];
	in[j-1] = in[i-1] - tempr;
	in[j  ] = in[i  ] - tempi;
	in[i-1] = in[i-1] + tempr;
	in[i  ] = in[i  ] + tempi;
      }
      tempr=wr;
      wr	+=tempr*wpr - wi   *wpi;
      wi	+=wi   *wpr + tempr*wpi;
    }
    mmax=istep;
  }
  return;
}

void realft(vector<double> &in)
{
  const int n = in.size()/2;		//the size of in must be double even
  int i1,i2,i3,i4;
  const double theta=	M_PI/n;
  const double wpr  =	-2*pow(sin(0.5*theta) , 2);
  const double wpi  =	sin(theta);
  double wr   = 1+wpr;
  double wi   = wpi;
  double h1i,h1r,h2i,h2r,wr_temp;
  double c1=0.5,c2=-0.5;
  fourier(in,1);

  for(int i=2;i<=n/2+1;++i)
  {
    i1=2*i-2;
    i2=i1+1;
    i3=2*n-i2+1;
    i4=i3+1;

    h1r= c1*(in[i1] + in[i3]);
    h1i= c1*(in[i2] - in[i4]);
    h2r=-c2*(in[i2] + in[i4]);
    h2i= c2*(in[i1] - in[i3]);

    in[i1] = h1r + wr*h2r - wi*h2i;
    in[i2] = h1i + wr*h2i + wi*h2r;
    in[i3] = h1r - wr*h2r + wi*h2i;
    in[i4] =-h1i + wr*h2i + wi*h2r;

    wr_temp=wr;
    wr	+=wr_temp*wpr - wi     *wpi;
    wi	+=wi     *wpr + wr_temp*wpi;
  }

  h1r  =in[0];
  in[0]=h1r+in[1];
  in[1]=h1r-in[0];
}

void sinft(vector<double> &in)
{
  if(in.size()%2!=0) 
  {
    cout<<"in.size  !=  out.size,error from sinft"<<endl;
    return;
  }
  const int n= in.size();		//the size of in must be even
  const double theta=	M_PI/n;
  const double wpr  =	-2*pow(sin(theta/2) , 2);
  const double wpi  =	sin(theta);
  double wr=1,wi=0;
  double wr_temp,y1,y2;
  in[0]=0;
  for(int i=1;i<=n/2;++i)
  {
    wr_temp=wr;
    wr	+=wr_temp*wpr - wi     *wpi;
    wi	+=wi     *wpr + wr_temp*wpi;
    y1	 =wi * (in[i] + in[n-i]);
    y2	 =0.5* (in[i] - in[n-i]);
    in[i]    =y1+y2;
    in[n-i]  =y1-y2;
  }
  realft(in);
  double sum=0;
  in[0]*=0.5;
  in[1] =0;
  for(int i=0;i<n;i+=2)
  {
    sum    +=in[i];
    in[i]  = in[i+1];
    in[i+1]= sum;
  }
  return;
}

