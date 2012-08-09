//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo....//
//                                                                    //
// Internal shock model for relativistic astrophysical jets           //
//                                                                    //
// Author: Omar Jamil                                                 //
// Contact: o.jamil@phys.soton.ac.uk                                  //
//          University of Southampton                                 //
//                                                                    //
//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo....//

#ifndef RADIATION_H
#define RADIATION_H
//Ref: Longair, '02 ~pg 250 - 262
#include <vector>
#include "shells.hh"
#include "container.hh"


//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
class Radiation
{
  
public:
  Radiation(); 
  //Radiation(const Radiation &c);
  ~Radiation();
  
  void synchSpectrum(std::vector<double> *, Container *);
  
  //Sunchrotron Power per unit volume per unit frequency
  double synchrotronEmiss(double &, double &, double &, double &);
  
  //Synchrotron absorption coefficient
  double synchrotronAbsor(double &, double &, double &, double &);

private:
  std::vector<double> *nu;
  std::vector<double> *iNu;
 
    
  
 
};
//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....

#endif // RADIATION_H
