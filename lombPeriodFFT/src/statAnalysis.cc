//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....//
// statistical analysis function such as Lomb-Scargle                  //
// Omar Jamil 2008                                                     //
// from routines used in Numerical recipes                             //
//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....//

#include "statAnalysis.hh"
#include <iostream>
#include <stdexcept>
#include <fstream>
#include <sstream>
#include <algorithm>

#define sqr(a) ((a)*(a))
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
StatAnalysis::StatAnalysis()
{}

//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
StatAnalysis::~StatAnalysis()
{}

//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
void StatAnalysis::getData(std::string &file, std::vector<double> &x, 
                           std::vector<double> &y)
  
{
  
  double timeData, fluxData, errorData;
  
  std::ifstream datafile(file.c_str());

  if (!datafile.good())
    throw std::runtime_error("Data file failed to load");
  
  while(!datafile.eof())
    {
      std::string line;
      getline(datafile, line);
      std::stringstream stream;
      
      if(line != "")
        {
          stream << line;
          stream >> timeData >> fluxData;
          
          x.push_back(timeData);
          y.push_back(fluxData);
        }
      
    }
  
  datafile.close();
  
}
//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
void StatAnalysis::writeData(std::string &file, std::vector<double> &x,
                             std::vector<double> &y)
{
  std::ofstream outfile(file.c_str());
  std::vector<double>::iterator i, j;
  
  for(i=x.begin(), j = y.begin(); 
      i!=x.end(); ++i, ++j)
    {
      outfile<<*i<<"\t"<<*j<<std::endl;
    }
  
}

//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
void StatAnalysis::writeData(std::string &file,
                             std::vector<double> &y)
{
  std::ofstream outfile(file.c_str());
  std::vector<double>::iterator i, j;
  
  for(j = y.begin(); 
      j!=y.end(); ++j)
    {
      outfile<<*j<<"\t";
      ++j;
      outfile<<*j<<std::endl;
    }
  
}
//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
void StatAnalysis::aveVar(std::vector<double> &data, double &ave, 
                          double &var)
{
  double s, ep;
  int j;
  
  int n = data.size();
  ave = 0.0;
  
  for(j=0; j<n; j++)
    {
      ave += data[j];
    }
  
  ave /= n;
  var = ep = 0.0;
  
  for(j=0; j<n; j++)
    {
      s = data[j] - ave;
      ep += s;
      var += s*s;
    }

  var=(var-ep*ep/n)/(n-1);

}

//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
void StatAnalysis::periodogram(std::vector<double> &x, 
                               std::vector<double> &y, 
                               const double ofac, const double hifac, 
                               std::vector<double> &px, 
                               std::vector<double> &py,
                               double &nout, int &jmax, double &prob)

//Given n data points with abscissas {x[0..n-1]} (need not be evenly spaced)
//and ordinates {y[0..n-1]}. Desired oversampling factor {ofac} (typical >= 4)
//the routine fills array {px[0..np-1]} with an increasing sequence of frequencies
//(not angular frequencies) up to {hifac} times the "average" Nyquist frequency and
//fills the array {py[0..np-1]} with the values of the Lomb normalized periodogram
//at those frequencies. The array {x} and {y} are not altered. {np}, the dimensions of
//{px} and {py}, must be at least as large as {nout}, the number of frequencies returned
//The routine also returns jmax such that {py[jmax]} is the maximum element on {py}, and 
//prob, an estimate of the significance of that maximum against the hypothesis of 
//random noise. A small value of {prob} indicates that a significant periodic signal
//is present.


{
  const double twoPi = 6.2831853071795586476;
  int i, j;
  double ave, c, cc, cwtau, effm, expy, pnow, pymax, s, ss, sumc, sumcy, sums, 
    sumsh, sumsy, swtau, var, wtau, xave, xdif, xmax, xmin, yy, arg, wtemp;
  
  int n = x.size();
  int np = px.size();
  
  
  std::vector<double> wi(n), wpi(n), wpr(n), wr(n);
  nout = 0.5 * ofac * hifac * n;
  
  if (nout > np) 
    {
      std::cerr<<"output array too short in period"<<std::endl;
    }
  
  aveVar(y, ave, var);
  
  if (var == 0.0) 
    {
      std::cerr<<"zero variance in period"<<std::endl;
    }

  xmax = xmin = x[0];
  
  for (j=0; j<n; j++)
    {
      if (x[j] > xmax) 
        {
          xmax=x[j];
        }
      
      if (x[j] < xmin)
        {
          xmin=x[j];
        }
            
    }
  
  xdif = xmax - xmin;
  xave = 0.5 * (xmax + xmin);
  
  pymax = 0.0;
  pnow = 1.0/(xdif*ofac);
  
  for (j=0; j<n; j++)
    {
      arg = twoPi*((x[j]-xave) * pnow);
      wpr[j] = -2.0 * sqr(sin(0.5 * arg));
      wpi[j] = sin(arg);
      wr[j] = cos(arg);
      wi[j] = wpi[j];
            
    }
  
  for (i=0; i<nout; i++)
    {
      px[i] = pnow;
      sumsh = sumc = 0.0;
      
      for (j=0; j<n; j++)
        {
          c = wr[j];
          s = wi[j];
          sumsh += s*c;
          sumc += (c-s) * (c+s);
        }
      
      wtau = 0.5 * atan2(2.0*sumsh, sumc);
      swtau = sin(wtau);
      cwtau = cos(wtau);
      sums = sumc = sumsy = sumcy = 0.0;
      
      for (j=0; j<n; j++)
        {
          s = wi[j];
          c = wr[j];
          ss = s * cwtau - c * swtau;
          cc = c * cwtau + s * swtau;
          sums += ss * ss;
          sumc += cc * cc;
          yy = y[j] - ave;
          sumsy += yy * ss;
          sumcy += yy * cc;
          wr[j] = ((wtemp = wr[j]) * wpr[j] - wi[j] * wpi[j]) + wr[j];
          wi[j] = (wi[j] * wpr[j] + wtemp * wpi[j]) + wi[j];
        }
      
      py[i] = 0.5 * (sumcy * sumcy/sumc + sumsy*sumsy/sums)/var;
      
      if(py[i] >= pymax) 
        {
          pymax = py[jmax=i];
        }
      
      pnow += 1.0/(ofac*xdif);
      
    }
  
  expy = exp(-pymax);
  effm = 2.0 * nout/ofac;
  prob = effm * expy;
  
  if (prob > 0.01) 
    {
      prob = 1.0 - pow(1.0-expy,effm);
    }

}

//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
void StatAnalysis::fasPer(std::vector<double> &x, 
                          std::vector<double> &y,
                          const double ofac, const double hifac,
                          std::vector<double> &wk1, 
                          std::vector<double> &wk2,
                          double &nout, int &jmax, double &prob
                          )
{
  const int Macc=4;
  int j, k, ndim, nfreq, nfreqt;
  double ave, ck, ckk, cterm, cwt, den, df, effm, expy, fac, fndim, hc2wt,
    hs2wt, hypo, pmax, sterm, swt, var, xdif, xmax, xmin;
  
  int n=x.size();
  int nwk=wk1.size();
  nout=0.5*ofac*hifac*n;
  nfreqt= int(ofac*hifac*n*Macc);
  nfreq=64;
  while (nfreq < nfreqt)
    {
      nfreq <<= 1;
    }
  ndim = nfreq << 1;
  
  if (ndim < nwk)
    {
      std::cerr<<"Workspace too small fasper"<<std::endl;
    }
  aveVar(y,ave,var);
  
  if (var == 0.0)
    {
      std::cerr<<"zero variance in fasper"<<std::endl;
    }
  
  
  xmin=x[0];
  xmax=xmin;
  for (j=1; j<n; j++)
    {
      if (x[j] < xmin)
        {
          xmin=x[j];
        }
      
      if (x[j] > xmax)
        {
          xmax=x[j];
        }
    }
  
  
  xdif=xmax-xmin;
  std::vector<double> wk1_t(ndim, 0.0);
  std::vector<double> wk2_t(ndim, 0.0);
  
  fac=ndim/(xdif*ofac);
  fndim=ndim;
  for (j=0;j<n;j++)
    {
      ck=fmod((x[j]-xmin)*fac,fndim);
      ckk=2.*(ck++);
      ckk=fmod(ckk, fndim);
      ++ckk;
      spread(y[j]-ave, wk1_t, ck, Macc);
      spread(1.0, wk2_t, ckk, Macc);
    }
  
  realFT(wk1_t, 1);
  realFT(wk2_t, 1);
  df=1./(xdif*ofac);
  pmax = -1.0;
  
  for (k=2,j=0;j<nout;j++,k+=2)
    {
      hypo = sqrt(wk2_t[k] * wk2_t[k] + wk2_t[k+1] * wk2_t[k+1]);
      hc2wt = 0.5 * wk2_t[k]/hypo;
      hs2wt = 0.5 * wk2_t[k+1]/hypo;
      cwt = sqrt(0.5+hc2wt);
      swt = SIGN(sqrt(0.5-hc2wt), hs2wt);
      den=0.5* n + hc2wt * wk2_t[k] + hs2wt * wk2_t[k+1];
      cterm = sqr(cwt * wk1_t[k] + swt * wk1_t[k+1])/den;
      sterm = sqr(cwt * wk1_t[k+1] - swt * wk1_t[k])/(n-den);
      wk1[j] = (j+1) * df;
      wk2[j] = (cterm + sterm)/(2.0 * var);
      
      if (wk2[j] > pmax)
        {
          pmax = wk2[jmax=j];
        }
    }
  
  expy = exp(-pmax);
  effm = 2.0 * nout/ofac;
  prob = effm*expy;
  
  if (prob > 0.01)
    {
      prob = 1.- pow(1.0-expy, effm);
    }
}

//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
void StatAnalysis::spread(const double y, std::vector<double> &yy, 
                          const double x, const int m)
{
  static int nfac[11]=
    {
      0,1,1,2,6,24,120,720,5040,40320,362880
    };
  
  int ihi, ilo, ix, j, nden;
  double fac;
  
  int n = yy.size();
  if(m > 10)
    {
      std::cerr<<"factorial too small in spread"<<std::endl;
    }
  
  ix=int(x);
  if (x == double(ix)) 
    {
      yy[ix -1] += y;
    }
  else
    {
      ilo = std::min(std::max(int(x-0.5*m),0),int(n-m));
      ihi=ilo+m;
      nden = nfac[m];
      fac = x-ilo-1;
      
      for (j=ilo+1; j<ihi; j++)
        {
          fac *= (x-j-1);
        }       
      yy[ihi-1] += y*fac / (nden*(x-ihi));
       
      for (j=ihi-1; j>ilo; j--)
        {
          nden = (nden/(j-ilo))*(j-ihi);
          yy[j-1] += y*fac/(nden*(x-j));
        }
    }
   
}
  
//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
void StatAnalysis::realFT(std::vector<double> &data, const int isign)
{
  int i, i1, i2, i3, i4;
  double c1 = 0.5, c2, h1r, h1i, h2r, h2i, wr, wi, wpr, wpi, wtemp, theta;
  
  int n = data.size();

  //check that n is power of 2
  bool f;
  unsigned int v = data.size()/2;
  f = (v & (v - 1)) == 0;
  if(!f)
    {
      double power = log10(n)/log10(2);
      power + 1;
      n = int(pow(2, int(power)));
    }
  
  theta = 3.141592653589793238/double(n>>1);
  if (isign == 1)
    {
      c2 = -0.5;
      four1(data, 1);
   
    }
  else
    {
      c2 = 0.5;
      theta = -theta;
    }
  
  wtemp = sin(0.5*theta);
  wpr = -2.0 * wtemp * wtemp;
  wpi = sin(theta);
  wr = 1.0 + wpr;
  wi=wpi;
  
  for (i=1; i <(n>>2); i++)
    {
  
      i2 = 1+(i1=i+i);
      i4 = 1+(i3=n-i1);
      h1r = c1 * (data[i1]+data[i3]);
      h1i = c1 * (data[i2]-data[i4]);
      h2r = -c2 * (data[i2]+data[i4]);
      h2i = c2 * (data[i1]-data[i3]);
      data[i1] = h1r+wr*h2r-wi*h2i;
      data[i2] = h1i+wr*h2i+wi*h2r;
      data[i3] = h1r-wr*h2r+wi*h2i;
      data[i4] = -h1i+wr*h2i+wi*h2r;
      wr = (wtemp=wr)*wpr-wi*wpi+wr;
      wi = wi*wpr+wtemp*wpi+wi;
    }

  if (isign == 1)
    {
      data[0] = (h1r=data[0])+data[1];
      data[1] = h1r-data[1];
    }
  else
    {
      data[0] = c1 * ((h1r=data[0])+data[1]);
      data[1] = c1*(h1r-data[1]);
      four1(data,-1);
    }
   
}
//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
void StatAnalysis::four1(std::vector<double> &data, const int isign)
{
  int n, mmax, m, j, istep, i;
  double wtemp,wr, wpr, wpi, wi, theta, tempr, tempi;

  int nn = data.size()/2;
    
  n = nn << 1;
  j=1;
  
  for (i=1; i<n; i+=2)
    {
      if (j > i)
        {
          std::swap(data[j-1], data[i-1]);
          std::swap(data[j], data[i]);
        }
      m = nn;
      while (m >= 2 && j >m)
        {
          j -= m;
          m >>= 1;
        }
      j +=  m;
    }
  
  mmax=2;
  while (n > mmax)
    {
      istep = mmax << 1;
      theta = isign*(6.28318530717959/mmax);
      wtemp = sin(0.5*theta);
      wpr = -2.0 * wtemp * wtemp;
      wpi = sin(theta);
      wr = 1.0;
      wi = 0.0;
      
      for (m=1; m<mmax; m+=2)
        {
          for (i=m; i<=n; i+=istep)
            {
              j = i+mmax;
              tempr = wr * data[j-1] - wi * data[j];
              tempi = wr * data[j] + wi * data[j-1];
              data[j-1] = data[i-1] - tempr;
              data[j] = data[i] - tempi;
              data[i-1] += tempr;
              data[i] += tempi;
                            
            }
          wr = (wtemp=wr)*wpr-wi*wpi+wr;
          wi = wi*wpr+wtemp*wpi+wi;
        }
      mmax = istep;
    }
  
}
//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
void StatAnalysis::fourierOnly(std::vector<double> &data, const int isign)
{
  int n, mmax, m, j, istep, i;
  double wtemp,wr, wpr, wpi, wi, theta, tempr, tempi;

  std::vector<double> data2;
  std::vector<double>::iterator d1, d2;
  
  for (d1 = data.begin();
       d1 != data.end(); ++d1)
    {
      data2.push_back(*d1);
      data2.push_back(0.0);
    }
  data.clear();
  
  for (d2 = data2.begin();
       d2 != data2.end(); ++d2)
    {
      data.push_back(*d2);
    }
  
  int nn = data.size()/2;
     
  n = nn << 1;
  j=1;
  
  for (i=1; i<n; i+=2)
    {
      if (j > i)
        {
          std::swap(data[j-1], data[i-1]);
          std::swap(data[j], data[i]);
        }
      m = nn;
      while (m >= 2 && j >m)
        {
          j -= m;
          m >>= 1;
        }
      j +=  m;
    }
  
  mmax=2;
  while (n > mmax)
    {
      istep = mmax << 1;
      theta = isign*(6.28318530717959/mmax);
      wtemp = sin(0.5*theta);
      wpr = -2.0 * wtemp * wtemp;
      wpi = sin(theta);
      wr = 1.0;
      wi = 0.0;
      
      for (m=1; m<mmax; m+=2)
        {
          for (i=m; i<=n; i+=istep)
            {
              j = i+mmax;
              tempr = wr * data[j-1] - wi * data[j];
              tempi = wr * data[j] + wi * data[j-1];
              data[j-1] = data[i-1] - tempr;
              data[j] = data[i] - tempi;
              data[i-1] += tempr;
              data[i] += tempi;
                            
            }
          wr = (wtemp=wr)*wpr-wi*wpi+wr;
          wi = wi*wpr+wtemp*wpi+wi;
        }
      mmax = istep;
    }
  
}
//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
Parameters::Parameters(std::string fn)
{
  try 
    {
      loadParametersFile(fn);
    } 
  catch (const std::exception &e) 
    {
      throw;
    }
    
}

//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
Parameters::~Parameters()
{}

//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
void Parameters::loadParametersFile(std::string filename)
  throw (const std::exception &)
{
  std::ifstream parameterFile(filename.c_str());
  if (!parameterFile.good())
    throw std::runtime_error("Parameter file failed to load");

  std::string label, dummy;
  std::cout<<"Reading the parameters file...."<<"\n";
  
  
  while(!parameterFile.eof())
    {
      parameterFile >> label;
      
      if(label == "input_file")
        {
          parameterFile >> dummy >> inputFile_;
        }
      else if(label == "output_file")
        {
          parameterFile >> dummy >> outputFile_;
        }
      else if(label == "analysis_choice")
        {
          parameterFile >> dummy >> analysisChoice_;
        }
      else if(label == "oversampling_factor")
        {
          parameterFile >> dummy >> overSampling_;
        }
      else if(label == "frequency_factor")
        {
          parameterFile >> dummy >> freqFactor_;
        }
                  
      
    }
  
}
//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
