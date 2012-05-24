//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo....//
//                                                                    //
// Internal shock model for relativistic astrophysical jets           //
//                                                                    //
// Author: Omar Jamil                                                 //
// Contact: o.jamil@phys.soton.ac.uk                                  //
//          University of Southampton                                 //
//                                                                    //
//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo....//

//putting the simulation together

#include <iostream>


#include "simulation.hh"
#include "evolution.hh"
#include "paraccess.hh"

//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
Simulation::Simulation()
{}

//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....

Simulation::~Simulation()
{
  delete time;
}
//Other classes destructors called before this one completes
//hence only the standard destructor here
//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....

void Simulation::run(int &duration, double &lum)
{
  try
    {
      bool useShFile = extUseShellFile;
      //deque containing all the shells
      //shells = new Container();
      //deque that will contain only the active shells
      liveShells = new Container();
      
      
      //initial population of the shells for randomly generated
      //shells. Written to shell file.
      if(!useShFile)
        {
          liveShells->populateWithShells(duration, lum);
        }
      
      //sorting the shells according to injection times
      //stable_sort(sh_dist->begin(),
      //            sh_dist->end(), TIfunc());
      
      //Evolution initialized with 'all' and 'active' deques
      time = new Evolution(liveShells);
            
      bool useIncreaseTStep = extIncreaseTimeStep;
      
  
      if(useIncreaseTStep)
        {
          time->evolveIncrSampling(duration);
        }
      else
        {
          time->evolve(duration);
        }
      
    }
  catch(std::bad_alloc)
    {
      std::cerr<<"Memory exhauseted in run;simulation.cc\n";
    }
  
  
}

//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....


