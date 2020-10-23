//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo....//
//                                                                    //
// Internal shock model for relativistic astrophysical jets           //
//                                                                    //
// Author: Omar Jamil                                                 //
// Contact: o.jamil@phys.soton.ac.uk                                  //
//          University of Southampton                                 //
//                                                                    //
//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo....//


#ifndef NUMERICAL_HH
#define NUMERICAL_HH

#include <vector>
#include <cmath>

//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
class Numerical
{
  
public:
  Numerical(); 	//the default constructor
  //Numerical(const Numerical &c);
  ~Numerical();
  
  //Pointer function; takes a const a double and returns a double
  typedef double (*pfn)(const double);
  //choosing which method to integrate with
  typedef double choose(pfn, const double, const double, const int);

  inline double rombergIntegration(pfn pf, const double &a, const double &b  
                            /*,choose */)
  {
    return romberg(pf,a, b/*,choose */);
  }
  

  inline double gammaFunction(const double &d)
  {
    return exp(gammalnFunction(d));
  }
  
  //protected:
  double random(int &);

private:
  
  //double choose(pfn, const double, const double, const int);
  
  //polynomial interpolation
  void polint(const std::vector<double> &, const std::vector<double> &, 
              const double , double &, double &);

  //integration function used in romberg; good for infinite limits 
  double midinf(pfn , const double , 
                const double , const int );
  
  double midexp(pfn , const double , 
                const double , const int );
  
  //the trapezium function for integration
  double trapz(pfn, const double , 
                const double , const int );
  
  //romberg method for integration
  //calls polint and choice of midinf or trapz
  double romberg(pfn ,const double , const double  
                 /*,choose */);
  
  double gammalnFunction(const double &);

  
    
  //change of limits for midinf (used for infinite limits)
  double funcInf(double (funk)(const double), const double x)
  {
    return funk(1.0/x)/(x*x);
  }
  
  double funcExp(double (funk)(const double), const double x)
  {
    return funk(-log(x))/x;
  }
};
//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....

#endif // NUMERICAL_HH
