//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo....//
//                                                                    //
// Internal shock model for relativistic astrophysical jets           //
//                                                                    //
// Author: Omar Jamil                                                 //
// Contact: o.jamil@phys.soton.ac.uk                                  //
//          University of Southampton                                 //
//                                                                    //
//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo....//


#ifndef CONTAINER_HH
#define CONTAINER_HH

#include <list>
#include <cmath>
#include <vector>

#include "shells.hxx"


//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
class Container : public std::list<Shells*>
{
  
public:
 
  //friend class Shells;
  
  // class to hold the collection of shells
  Container();	//the default constructor
  //Container(const Container &c) {}
  ~Container();

  //populating the Container with shells
  void populateWithShells(int &, double &);
 
  //update the shells dimension for the elapsed time 
  void updateShellParams(double &);
  
  void coutContainerElements();
  
  void majorEjection(double &, double &);
  
  //routine to prevent the shells overlapping or overtaking
  void overlapPrevention(double &);
  
  //This reads the shell data file and then writes the shell
  //properties to the vectors.
  void shellData(std::vector<double> &, std::vector<double> &,
                 std::vector<double> &, std::vector<double> &);
  
  //"Random" sampling from gaussian dist of given mean and std
  double gaussian(const double &, const double &, double, double);
  
  bool fileExists(const std::string& );
  
  void injectShell(double &, double &, double &, double &);
  
  void dumpActiveShells();
  
  void readActiveShells();
  
  //function to read the lightcurve data (currently redundant)
 //  void getData();
  
//   //fucntion to read shell data (produced from a lightcurve)
//   //(currently redundant)
//   void shellData();
 
protected:
  void outOfStore()
  {
    std::cerr<<"Operator new failed in Container.cc \n";
    throw std::bad_alloc();
  }
  
private:
  
  //vectors to store the lightcurve data
  // std::vector<double> *timeAxis;
//   std::vector<double> *fluxAxis;
//   //vectors to store Injection time, Mass,  shell Gamma (BLF), length
//   std::vector<double> *shITime;
//   std::vector<double> *shM;
//   std::vector<double> *shGam;
//   std::vector<double> *shLength;
     
  
};
//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....









#endif // CONTAINER_HH
