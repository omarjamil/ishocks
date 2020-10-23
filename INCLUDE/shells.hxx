//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo....//
//                                                                    //
// Internal shock model for relativistic astrophysical jets           //
//                                                                    //
// Author: Omar Jamil                                                 //
// Contact: o.jamil@phys.soton.ac.uk                                  //
//          University of Southampton                                 //
//                                                                    //
//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo....//

#ifndef SHELLS_HH
#define SHELLS_HH

#include <cmath>
#include <vector>
#include <iostream>

//#include "parameters.h"

//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....


class Shells
{

public:

  //  constructor
  Shells(double &, double &, double &, double &, double &);

  //Default generated copy constructor in use
  //Shells(const Shells&); //copy constructor

  ~Shells();

  //the initial shell distribution


  void shellInitialization();

  //setting the shell dimensions x_lower, x_upper, vol
  //at various time steps
  void setDimensions(double &);

  //setting the volume and area for the merged shell
  void setInitialMergedVA();

  //setting the lower and upper x coordinates
  //separate from setDimensions as this one is
  //used by merger formed shells
  void setLowerAndUpperX(double &, double &);

  void setMergeCount(int &);

  //the shell expansion velocity
  void setExpansionBeta();

  //powerlaw Normalization calculated with Eint in units of GeV (Longair)
  void powerlawNorm();

  //electron powerlaw distribution
  void powerLaw();

  //Doppler effect factor
  double dopplerFactor();

  //Changing the powerlaw normalization based on the change in the shell
  //volume. This should account for the adiabatic losses.
  void adiabaticLosses(double &);

  //setting the Magentic energy density
  void initialBEnergyDensity();

  //setting magnetic field parameters after adiabatic losses
  void setMagParams(double &);

  //if the shells were expanding into vacuum i.e. no work done.
  void inVacuumExpansion();

  //for slower energization of merged shells
  void slowEnergization(double &);

  void shockZone(double &);

  //radiative losses (experimental)
  void radiativeLosses(double &);

  //inline functions, does not necessarily make it
  //quicker but more of a suggestion for the compiler

  //get the time of injection
  inline double getTimeOfInjection()
  {
    return ti;
  }

  //if a shell then will return true i.e. 1
  inline bool getShellId()
  {
    return shell;
  }

  inline void setShellId(bool &id)
  {
    shell = id;
  }

  inline double getShellGamma()
  {
    return gamma;
  }

  inline double getShellMass()
  {
    return mass;
  }

  inline double getShellVolume()
  {
    return vol;
  }

  inline void setShellVolume(double &V)
  {
    vol = V;
  }

  //returns the curved surface area
  inline double getShellArea()
  {
    return area;
  }

  inline void setShellArea(double &A)
  {
    area = A;
  }

  inline double getShellWidth()
  {

    return (x_u - x_l);
  }

  //distance from the source to outer edge of the shell
  inline double getOuterRadius()
  {
    return x_u;
  }

  inline void setOuterRadius(double &ouRad)
  {
    x_u = ouRad;
  }

  //distance from the source to inner edge of the shell
  inline double getInnerRadius()
  {
    return x_l;
  }

  inline void setInnerRadius(double &innRad)
  {
    x_l = innRad;
  }

  //Radius of the inner edge of the shell i.e. the lower
  //face of the frustum radius.
  inline double getShellInRad()
  {
    return r_l;
  }

  inline void setShellInRad(double &Rl)
  {
    r_l = Rl;
  }

  inline double getShellOuRad()
  {
    return r_u;
  }

  inline void setShellOuRad(double &Ru)
  {
    r_u = Ru;
  }

  inline double getShellBeta()
  {
    return sqrt(1. - (1. / pow(gamma, 2)));
  }

  inline void setInternalEnergy(double &intEner)
  {
    if(intEner <= 0.0)
      {
        ienergy = 0.0;
      }
    else
      {
        ienergy = intEner;
      }
  }


  inline double getInternalEnergy()
  {
    return ienergy;
  }

  inline double getThermalEnergy()
  {
    return shellThermEner;
  }

  inline void setThermalEnergy(double &Eth)
  {
    shellThermEner = Eth;
  }


  inline double getLocation()
  {
    // double loc = x_l + (getShellWidth()/2.);

//     return loc;
    return shellCentre;
  }

  inline void setLocation(double &loc)
  {
    shellCentre = loc;

  }


  inline void setIdFalse()
  {
    shell = false;
  }

  inline int getMergeCount()
  {
    return mergeCount;
  }


  inline double getExpansionBeta()
  {
    return betaExpansion;
  }

  inline void setExpansionBeta(double &betE)
  {
    betaExpansion = betE;
  }

   inline double getPLawNorm()
  {

    return pNorm;

  }

  inline void setPLawNorm(double &p)
  {
    pNorm = p;
  }

  inline void nuWrite(double &iNu)
  {
    iNuGrid->push_back(iNu);
  }

  inline void copyINu(std::vector<double> &iN)
  {
    (*iNuGrid) = iN;
  }

  inline void copyTau(std::vector<double> &ta)
  {
    (*tau) = ta;
  }

  inline void nuClear()
  {
    iNuGrid->clear();
  }

  inline void tauWrite(double &tauVal)
  {
    tau->push_back(tauVal);
  }

  inline void tauClear()
  {
    tau->clear();
  }

  inline std::vector<double> iNuVals()
  {
    return (*iNuGrid);

  }

  inline std::vector<double> tauVals()
  {
    return (*tau);
  }

  inline std::vector<double> eleDistVals()
  {
    return (*gammaGrid);
  }

  inline std::vector<double> gammaVals()
  {
    return (*egamma);
  }


  inline double getGammaMin()
  {
    return egammaMin;
  }

  inline void setGammaMin(double &eGm)
  {
    egammaMin = eGm;
  }

  inline double getGammaMax()
  {
    return egammaMax;
  }

  inline void setGammaMax(double &gMax)
  {
    egammaMax = gMax;
  }

  inline double getBEneDens()
  {
    return BEDens;
  }

  inline void setBEneDens(double &bED)
  {
    BEDens = bED;
  }

  inline double getMagPressure()
  {
    return magPressure;

  }

  inline void setMagneticPressure(double &magEneDens)
  {
    magPressure = 1./3. * magEneDens;
  }

  inline void setPrevTime(double &prevT)
  {
    tPrev = prevT;
  }

  inline double getPrevTime()
  {
    return tPrev;
  }


  //The for reaching maximum internal energy
  inline double getTEnergize()
  {
    return tEnergize;
  }

  inline void setTEnergize(double &tE)
  {
    tEnergize = tE;
  }

  inline double getIntEfrac()
  {
    return fracEInt;
  }

  inline void setIntEfrac(double &Efrac)
  {
    fracEInt = Efrac;
  }

  inline double getIntEavail()
  {
    return intEavail;
  }

  inline void setIntEavail(double &inEAvail)
  {
    intEavail = inEAvail;
  }

  inline void shockShell()
  {
    shocked = 1;
  }

  inline bool shellShocked()
  {
    return shocked;
  }

  //used when starting from previous simulation shells.
  inline void setShellShocked(bool &sss)
  {
    shocked = sss;
  }


private:

  double ti;
  double mass;
  double gamma;
  double ienergy;
  double width; //shell length along the jet
  double betaExpansion;
  double x_l; //Distance from source to the inner part of the shell
  double x_u;
  double r_l; //Shell's inner radius
  double r_u;
  double vol;
  double area;
  double egammaMin;
  double egammaMax;
  int mergeCount;
  bool shell;
  double pNorm;
  double BEDens;
  double magPressure;
  double shellCentre;
  //the spectrum container
  std::vector<double> *iNuGrid;
  std::vector<double> *gammaGrid;
  std::vector<double> *egamma;

  std::vector<double> *tau;
  double tPrev;
  double shellThermEner;
  double tEnergize; //time for reaching max internal energy
  double fracEInt;
  double intEavail;
  bool shocked;

};
//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....


#endif // SHELLS_HH
