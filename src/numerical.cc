//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo....//
//                                                                    //
// Internal shock model for relativistic astrophysical jets           //
//                                                                    //
// Author: Omar Jamil                                                 //
// Contact: o.jamil@phys.soton.ac.uk                                  //
//          University of Southampton                                 //
//                                                                    //
//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo....//

#include <cmath>
#include <iostream>
#include <vector>


#include "numerical.hh"


//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
Numerical::Numerical()
{}

//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....

Numerical::~Numerical()
{}

//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
void Numerical::polint(const std::vector<double> &xa, 
                       const std::vector<double> &ya, 
                       const double x, double &y, double &dy)
{
  int i, m, ns=0;
  double den, dif, dift, ho, hp, w;
  
  int n=xa.size();
  std::vector<double> c(n), d(n);
  dif = fabs(x-xa[0]);
  
  for (i = 0; i < n; i++)
    {
      if((dift=fabs(x-xa[i])) < dif )
        {
          ns =i;
          dif=dift;
        }
      c[i]=ya[i];
      d[i]=ya[i];
    }
  
  y=ya[--ns];
  
  for (m = 1; m < n; m++)
    {
      for (i = 0; i < n-m; i++)
        {
          ho=xa[i]-x;
          hp=xa[i+m]-x;
          w=c[i+1]-d[i];
          if ((den=ho-hp) == 0.0) {std::cout<<"Error in routine polint"<<"\n";}
          den=w/den;
          d[i]=hp*den;
          c[i]=ho*den;
          
        }
      y += (dy=(2*(ns+1) < (n-m) ? c[ns+1] : d[ns--]));
      
    }
  
}
//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
double Numerical::trapz(pfn func, const double a, 
                        const double b, const int n)
{
  double x, tnm, sum, del;
  double s=0.0;
  int it, j;
  
  if(n==1){
    return (s =0.5*(b-a)*(func(a)+func(b)));
  }
  else {
    for (it=1, j=1; j<n-1; j++) it <<=1;
    tnm=it;
    del=(b-a)/tnm;
    x=a+0.5*del;
    for(sum=0.0,j=0;j<it;j++,x+=del)sum += func(x);
    s=0.5*(s+(b-a)*sum/tnm);
    return s;
  }
  
  
}

//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....

double Numerical::midinf(pfn funk, const double aa, 
                         const double bb, const int n)
{

  double x, tnm, sum, del, ddel, b, a;
  double s=0.0;
  int it, j;
  
  b=1./aa;
  a=1./bb;
  
  if(n ==1)
    {
      return (s=(b-a)*funcInf(funk, 0.5*(a+b)));
      
    }
  else 
    {
      for (it=1, j=1; j<n-1;j++) it *=3;
      tnm=it;
      del=(b-a)/(3.0*tnm);
      ddel = del+del;
      x=a+0.5*del;
      sum=0.0;
      
      for (j=0;j<it;j++)
        {
          sum += funcInf(funk, x);
          x += ddel;
          sum += funcInf(funk, x);
          x += del;
          
        }
      return (s=(s+(b-a)*sum/tnm)/3.0);
      
    }
   

}
//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
double Numerical::midexp(pfn funk, const double aa, 
                         const double bb, const int n)
{

  double x, tnm, sum, del, ddel, b, a;
  double s=0.0;
  int it, j;
  
  b=exp(-aa);
  a=1.;
  
  if(n ==1)
    {
      return (s=(b-a)*funcExp(funk, 0.5*(a+b)));
      
    }
  else 
    {
      for (it=1, j=1; j<n-1;j++) it *=3;
      tnm=it;
      del=(b-a)/(3.0*tnm);
      ddel = del+del;
      x=a+0.5*del;
      sum=0.0;
      
      for (j=0;j<it;j++)
        {
          sum += funcExp(funk, x);
          x += ddel;
          sum += funcExp(funk, x);
          x += del;
          
        }
      return (s=(s+(b-a)*sum/tnm)/3.0);
      
    }
   

}

//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....

double Numerical::romberg(pfn func, const double a, const double b 
                          /*,choose workFunc*/)
{
  const int JMAX=20, JMAXP=JMAX+1, K=5;
  const double EPS=3.e-10;
  double ss, dss;
  std::vector<double> s(JMAX), h(JMAXP), s_t(K), h_t(K);
  int i, j;
  
  h[0] = 1.0;
  for (j = 1; j <= JMAX; j++)
    {
      //s[j-1]=workFunc(func, a, b, j);
      //s[j-1]=midinf(func, a, b, j);
      s[j-1]=trapz(func, a, b, j);
      //s[j-1]=midexp(func, a, b, j);
      
      if (j >= K)
        {
          for (i = 0; i < K; i++)
            {
              h_t[i]=h[j-K+i];
              s_t[i]=s[j-K+i];
              
            }
          polint(h_t, s_t, 0.0, ss, dss);
          //std::cout<<fabs(dss)<<"\t"<<fabs(ss)<<"\n";
          if (fabs(dss) <= EPS*fabs(ss)) {return ss;}
          
          
          
        }
      h[j]=h[j-1]/10.0; //use with midinf
      //h[j]=0.25*h[j-1]; //use with trapz
      
    }
  std::cout<<"Too many steps in routine qromb"<<"\n";
  return 0.0;
  
}
//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
double Numerical::gammalnFunction(const double &xx)
{
  int j;
  double x, y, tmp, ser;
  const double cof[6]={76.18009172947146, -86.50532032941677,
                              24.01409824083091, -1.231739572450155,
                              0.1208650973866179e-2, -0.5395239384953e-5};
  
  y=x=xx;
  tmp=x+5.5;
  tmp -= (x+0.5)*log(tmp);
  ser=1.000000000190015;
  
  for(j=0;j<6;j++)
    {
      ser += cof[j]/++y;
    }
  
  return -tmp+log(2.50662827746310005*ser/x);
  
  
}

//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
//From Numerical Recipes ran2
double Numerical::random(int &idum)
{
  const int im1=2147483563, im2=2147483399;
  const int ia1=40014, ia2=40692, iq1=53668, iq2=52774, ir1=12211;
  const int ir2=3791, ntab=32, imm1=(im1-1);
  const int ndiv=(1+imm1/ntab);
  const double eps=1.2e-7, RNMX=(1.0-eps), am=(1./double(im1));
  int j, k;
  static int idum2=123456789;
  static int iy=0;
  static std::vector<int> iv(ntab);
  double temp;
  
  if (idum <= 0) {
    if (-(idum) < 1) idum=1;
    else idum = -(idum);
    idum2=(idum);
    for (j=ntab+7;j>=0;--j) {
      k=(idum)/iq1;
      idum=ia1*(idum-k*iq1)-k*ir1;
      if (idum < 0) idum += im1;
      if (j < ntab) iv[j] = idum;
    }
    iy=iv[0];
  }
  k=(idum)/iq1;
  idum=ia1*(idum-k*iq1)-k*ir1;
  if (idum < 0) idum += im1;
  k=idum2/iq2;
  idum2=ia2*(idum2-k*iq2)-k*ir2;
  if (idum2 < 0) idum2 += im2;
  j=iy/ndiv;
  iy=iv[j]-idum2;
  iv[j] = idum;
  if (iy < 1) iy += imm1;
  if ((temp=am*iy) > RNMX) return RNMX;
  else return temp;
  
}

//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
