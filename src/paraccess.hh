//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo....//
//                                                                    //
// Internal shock model for relativistic astrophysical jets           //
//                                                                    //
// Author: Omar Jamil                                                 //
// Contact: o.jamil@phys.soton.ac.uk                                  //
//          University of Southampton                                 //
//                                                                    //
//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo....//

//The following file contains the collection of parameters
//for the simulation. Including this file should make these available.

#ifndef PARACCESS_HH
#define PARACCESS_HH

extern bool extLoadPrevSim;
extern bool extDumpActive;
extern double extJetLum;
extern double extJetOpenAngle;
extern double extJetViewAngle;
extern double extSourceDistance;
extern double extEjTimeGap;
extern int extSimDuration;
extern double extGamMax;
extern double extGamMin;
extern double extShWidth;
extern double extShThermalEne;
extern double extShIntKine;
extern double extShIntMag;
extern double extPowLawInd;
extern double exteGammaMin;
extern double exteGammaMax;
extern double extnuMin;
extern double extnuMax;
extern int extnuPoints;
extern bool extindivFreq;
extern double extnu1;
extern double extnu2;
extern double extSynchConstAp;
extern double extSynchConstBp;
extern bool extUseShellFile;
extern std::string extShellFile;
extern bool extUseDataFile;
extern std::string extDatFile;
extern std::string extResultsFile;
extern int extEveryTStep;
extern bool extIncreaseTimeStep;
extern double extStepResolution;
extern double extTotalDuration;
extern bool extInVacuum;
extern bool extInjIntEne;
extern double extRelMassFrac;
extern bool extSlowEnerg;
extern bool extUseResFile;
extern std::string extLCFile;
extern bool extFinalTStep;
extern double extShockLoc;
//****TEMPORARY??****
extern double extUp1;
extern double extUp2;
extern double extUp3;
extern double extUp4;
extern double extUp5;
extern bool extSplitLC;
//****??TEMPORARY****
#endif // PARACCESS_HH
