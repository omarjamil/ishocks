//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo....//
//                                                                    //
// Internal shock model for relativistic astrophysical jets           //
//                                                                    //
// Author: Omar Jamil                                                 //
// Contact: o.jamil@phys.soton.ac.uk                                  //
//          University of Southampton                                 //
//                                                                    //
//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo....//

//the run file containing the main function.

#include <iostream>
#include <iterator>
#include <iostream>
#include <ctime>

#include "simulation.hxx"
#include "parameters.hxx"
#include "paraccess.hxx"


int main()
{
  
  
  // int duration = 200;
//   double luminosity = 2e24;
  
  std::time_t timeStart, timeEnd;
  time(&timeStart);
  
  Simulation sim;
  
  //simDuration and jetLum accessed from paraccess.h
  //that reads the parameter file in.
   
  sim.run(extSimDuration, extJetLum);
  
  time(&timeEnd);
  double diff = difftime(timeEnd, timeStart) / 60.;
  std::cout<<"CPU time usage for the Simualtion [minutes]: "<<diff<<std::endl;
    
  return 0;
  
}

