//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo....//
//                                                                    //
// Internal shock model for relativistic astrophysical jets           //
//                                                                    //
// Author: Omar Jamil                                                 //
// Contact: o.jamil@phys.soton.ac.uk                                  //
//          University of Southampton                                 //
//                                                                    //
//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo....//

// The following file links the parameters and makes them available
// to all other classes. The variables declared here need to declared
// with extern keyword in other classes. The following are the
// definitions.

#include "parameters.hxx"

Parameters params("parameters.par");
bool extLoadPrevSim = params.loadPrevSim();
bool extDumpActive = params.dumpActive();
double extJetLum = params.jetLuminosity();
double extJetOpenAngle = params.openingAngle();
double extJetViewAngle = params.viewingAngle();
double extSourceDistance = params.sourceDistance();
double extEjTimeGap = params.ejectionGap();
int extSimDuration = params.duration();
double extGamMax = params.GammaMax();
double extGamMin = params.GammaMin();
double extShWidth = params.shellWidth();
double extShThermalEne = params.shellThermal();
double extShIntKine = params.shellIntKine();
double extShIntMag = params.shellMagnetic();
double extPowLawInd = params.ePowerLawIndex();
double exteGammaMin = params.eGammaMin();
double exteGammaMax = params.eGammaMax();
double extnuMin = params.nuMin();
double extnuMax = params.nuMax();
int extnuPoints = params.nuPoints();
int extgammaPoints = params.gammaPoints();
bool extindivFreq = params.indivFreq();
double extnu1 = params.nu_1();
double extnu2 = params.nu_2();
bool extUseShellFile = params.useShellFile();
double extSynchConstAp = params.synchConstA_P();
double extSynchConstBp = params.synchConstB_P();
std::string extShellFile = params.shellFileName();
bool extUseDataFile = params.useDataFile();
std::string extDatFile = params.dataFileName();
std::string extResultsFile = params.resultsFileName();
int extEveryTStep = params.everyTStep();
bool extIncreaseTimeStep = params.increaseTimeStep();
double extStepResolution = params.stepResolution();
double extTotalDuration = params.totalDuration();
bool extInVacuum = params.inVacuum();
bool extInjIntEne = params.injIntEne();
double extRelMassFrac = params.relMassFrac();
bool extSlowEnerg = params.slowEnergization();
bool extUseResFile = params.useResultsFile();
std::string extLCFile = params.lcFile();
bool extFinalTStep = params.finalTStep();
double extShockLoc = params.shocLoc();
//****TEMPORARY??****
double extUp1 = params.Up1();
double extUp2 = params.Up2();
double extUp3 = params.Up3();
double extUp4 = params.Up4();
double extUp5 = params.Up5();
bool extSplitLC = params.SplitLC();
//****??TEMPORARY****
