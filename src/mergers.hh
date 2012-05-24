//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo....//
//                                                                    //
// Internal shock model for relativistic astrophysical jets           //
//                                                                    //
// Author: Omar Jamil                                                 //
// Contact: o.jamil@phys.soton.ac.uk                                  //
//          University of Southampton                                 //
//                                                                    //
//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo....//

#ifndef MERGERS_HH
#define MERGERS_HH

#include "shells.hh"
#include "container.hh"

//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....

//shells merger class

class Mergers 
{
  
public:
  
  //the default constructor
  Mergers(Shells *, Shells *);
  //  Mergers(const Mergers &c) {}
  ~Mergers();
  
  void doMerger(Container *, double &);
  
  double inNewShellWidth(double &, double &, 
                         double &, double &, 
                         double &);
  
  double ouNewShellWidth(double &, double &, 
                         double &, double &, 
                         double &);
  
  double mergedShellIntEner(double &, double &,
                            double &, double &,
                            double &, double &,
                            double &);
  
  
  double getInnerShell()
  {
    return innerShell->getTimeOfInjection();
  }
  
  double getOuterShell()
  {
    return outerShell->getTimeOfInjection();
  }
  
  void mergeLocCount(double &);
  
  //protected:
  
private:
  Shells *innerShell;
  Shells *outerShell;
  double i;
  double m;
  double g;
  bool sh;
  
  
};
//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....



#endif // MERGERS_HH
