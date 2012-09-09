//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo....//
//                                                                    //
// Internal shock model for relativistic astrophysical jets           //
//                                                                    //
// Author: Omar Jamil                                                 //
// Contact: o.jamil@phys.soton.ac.uk                                  //
//          University of Southampton                                 //
//                                                                    //
//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo....//

//putting the simulation together.

#ifndef SIMULATION_HH
#define SIMULATION_HH


#include <algorithm>
#include <functional>

#include "funcobj.hxx"
#include "container.hxx"
#include "evolution.hxx"

//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
class Simulation
{
  
public:
  Simulation(); 	//the default constructor
  //Simulation(const Simulation &c);
  ~Simulation();
  
  void run(int &, double &);
  
  //protected:
  
  
  
private:
  
  //Container *shells;
  Container *liveShells;
  Evolution *time;
   
    
};

//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....

#endif // SIMULATION_HH
