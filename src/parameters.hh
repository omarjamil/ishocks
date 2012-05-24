//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo....//
//                                                                    //
// Internal shock model for relativistic astrophysical jets           //
//                                                                    //
// Author: Omar Jamil                                                 //
// Contact: o.jamil@phys.soton.ac.uk                                  //
//          University of Southampton                                 //
//                                                                    //
//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo....//

#ifndef PARAMETERS_HH
#define PARAMETERS_HH

#include <iostream>
#include <fstream>
#include <stdexcept>
#include <string>



//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
class Parameters
{
  
public:
  Parameters(std::string); 	//the default constructor
  //Parameters(const Parameters &c);
  ~Parameters();
  
  void loadParametersFile(std::string) throw (const std::exception &);

  bool loadPrevSim()
  {
    return loadPrevSim_;
  }
  
  bool dumpActive()
  {
    return dumpActive_;
  }
   
  double jetLuminosity()
  {
    return jetLuminosity_;
  }
  
  double openingAngle()
  {
    return openingAngle_;
  }

  double viewingAngle()
  {
    return viewingAngle_;
  }
  
  double sourceDistance()
  {
    return (sourceDistance_ * 3.08568025e16);
  }
   
  double ejectionGap()
  {
    return ejectionGap_;
  }
  
 
  int duration()
  {
    return duration_;
  }
  
  double GammaMax()
  {
    return GammaMax_;
  }

  double GammaMin()
  {
    return GammaMin_;
  }
  
  double shellWidth()
  {
    return width_;
  }
  
  double shellThermal()
  {
    return thermalEne_;
  }
  
  double shellIntKine()
  {
    return EintKine_;
  }
  
  double shellMagnetic()
  {
    return EMagnetic_;
  }
  
  double ePowerLawIndex()
  {
    return pLawIndex_;
  }

  double eGammaMin()
  {
    return gMin_;
  }
  
  double eGammaMax()
  {
    return gMax_;
  }
  
  double nuMin()
  {
    return nuMin_;
  }
  
  double nuMax()
  {
    return nuMax_;
  }
  
  int nuPoints()
  {
    return nuPoints_;
  }
  
  bool indivFreq()
  {
    return individual_frequencies_;
  }
  
  double nu_1()
  {
    return nu_1_;
  }
  
  double nu_2()
  {
    return nu_2_;
  }
  
  
  double synchConstA_P()
  {
    return synchConstap();
    
  }
  
  double synchConstB_P()
  {
    return synchConstbp();
    
  }

  bool useShellFile()
  {
    return useShellFile_;
  }
  
  std::string shellFileName()
  {
    return shellFile_;
  }
  
  bool useDataFile()
  {
    return useDataFile_;
  }

  std::string dataFileName()
  {
    return dataFile_;
  }
  
  std::string resultsFileName()
  {
    return resultsFile_;
  }

  int everyTStep()
  {
    return everyTStep_;
  }
  
  bool increaseTimeStep()
  {
    return increaseTimeStep_;
  }

  double stepResolution()
  {
    return stepResolution_;
  }
  
  double totalDuration()
  {
    return totalDuration_;
  }
  
  bool inVacuum()
  {
    return inVacuum_;
  }
  
  bool injIntEne()
  {
    return injIntEne_;
  }

  double relMassFrac()
  {
    return relMassFrac_;
    
  }
  
  bool slowEnergization()
  {
    return slowEnergization_;
  }
  
  bool useResultsFile()
  {
    return useResultsFile_;
  }
  
  std::string lcFile()
  {
    return lcFile_;
  }

  bool finalTStep()
  {
    return finalTStep_;
  }
  
  double shocLoc()
  {
    return shocLoc_;
  }
  //****TEMPORARY??****
  double Up1()
  {
    return up1;
  }

  double Up2()
  {
    return up2;
  }
  
  double Up3()
  {
    return up3;
  }
  
  double Up4()
  {
    return up4;
  }
  
  double Up5()
  {
    return up5;
  }

  bool SplitLC()
  {
    return split;
  }
  //****??TEMPORARY****
private:
  bool loadPrevSim_;
  bool dumpActive_;
  double jetLuminosity_;
  double openingAngle_;
  double viewingAngle_;
  double sourceDistance_;
  double ejectionGap_;
  int duration_;
  double GammaMax_;
  double GammaMin_;
  double width_;
  double thermalEne_;
  double EintKine_;
  double EMagnetic_;
  double pLawIndex_;
  double gMin_;
  double gMax_;
  double nuMin_;
  double nuMax_;
  int nuPoints_;
  bool individual_frequencies_;
  double nu_1_;
  double nu_2_;
  bool useShellFile_;
  bool increaseTimeStep_;
  bool useDataFile_;
  bool useResultsFile_;
  std::string shellFile_;
  std::string dataFile_;
  std::string resultsFile_;
  int everyTStep_;
  std::string lcFile_;
  double stepResolution_;
  double totalDuration_;
  bool inVacuum_;
  bool injIntEne_;
  double relMassFrac_;
  bool slowEnergization_;
  double synchConstap();
  double synchConstbp();
  bool finalTStep_;
  double shocLoc_;
 //****TEMPORARY??****
  bool split;
  double up1;
  double up2;
  double up3;
  double up4;
  double up5;
  //****??TEMPORARY****
  
};
//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....

#endif // PARAMETERS_HH
