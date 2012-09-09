//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo....//
//                                                                    //
// Internal shock model for relativistic astrophysical jets           //
//                                                                    //
// Author: Omar Jamil                                                 //
// Contact: o.jamil@phys.soton.ac.uk                                  //
//          University of Southampton                                 //
//                                                                    //
//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo....//

// Time evolution of the shells

#ifndef EVOLUTION_HH
#define EVOLUTION_HH


#include <vector>

#include "container.hxx"
#include "mergers.hxx"
//#include "parainclude.h"

//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....


class Evolution
{
public:
  Evolution(Container *);
  //Evolution(const Evolution &c);
  ~Evolution();
  //protected:

  //this runs the simulation based on time steps that are simply
  //"events" i.e. injection or merger.
  void evolve(int &);

  //this runs the simulation with increases samplings
  //appropriate for sampling shells without the need of
  //"event".
  void evolveIncrSampling(int &);


  void continueEvol(double &, Container *);


  //the min of all the collision times for active shells
  double minCollisionTime(Container *, std::vector<Mergers*> &);

  //calculating the collision time for two shells
  double collisionTime(Shells *, Shells *);

  //time til next injection
  double dtNextInjection(Shells *, Shells *);

  void results(Container *, std::vector<double> *, double &, int &);

  void results(Container *, std::vector<double> *, double &,
               std::vector<double> *, std::vector<double> *);

  //function to populate the frequencies container
  //logarithmic ranging
  void freqRange(std::vector<double> *);
 
  void avgLC(std::vector<double> *, std::vector<double> *);

  void finalTStep(Container *, std::vector<double> *,
                 double &);

  void finalTStepTau(Container *, std::vector<double> *,
                 double &);

  void tauWrite(Container *, std::vector<double> *);

  void splitLightcurve(Container *, double &);

private:
  //evolving the following Shells containers
  //being passed as pointers to Evolution class
  Container *activeShells;
  Mergers *merger;
  std::vector<double> *nu, *egamma, *comptonVect;

  //std::vector<double> *nuCont;
  //std::vector<double> *iNu;
  std::vector<double> *fluxNu1;
  std::vector<double> *fluxNu2;
  //std::vector<double> *fluxNu1Cont;
  //std::vector<double> *fluxNu2Cont;
  bool fileExists(const std::string& );
  

};
//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....

#endif // EVOLUTION_HH
