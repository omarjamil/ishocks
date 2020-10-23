//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo....//
//                                                                    //
// Internal shock model for relativistic astrophysical jets           //
//                                                                    //
// Author: Omar Jamil                                                 //
// Contact: o.jamil@phys.soton.ac.uk                                  //
//          University of Southampton                                 //
//                                                                    //
//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo....//

#include <cmath>
#include <iostream>
#include <sstream>
#include <string>
#include <iomanip>
#include <cstdlib>  // For srand() and rand()
#include <iterator>
#include <fstream>
#include <sstream>
#include <stdexcept>
//#include <memory>

#include "container.hxx"
#include "shells.hxx"
#include "physcon.hxx"
#include "paraccess.hxx"
#include "numerical.hxx"

#define sqr(a) ((a)*(a))

//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....

//shellsdeque class constructor
Container::Container():std::list<Shells*>()
{}

//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....

//shellsdeque class destructor
Container::~Container()
{

 

  for (iterator slist = begin();
       slist != end(); ++slist)
    {
      delete *slist;
    }
  
  //used in with getData() and shellData()
  //delete timeAxis;
  //delete fluxAxis;
  //delete shITime;
  //delete shM;
  //delete shGam;
  //delete shLength;
  
}


//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....

void Container::populateWithShells(int &tw, double &lum)
{
  //the time between injections
  double tv = extEjTimeGap;
  //the number of shells to injected
  // tw = total run-time of the simulation
  double n = round(tw/tv);
  //int n=10; //max number of shells in one time bin
  //injection time of the shell
  double ti=0.0;
    
  //Gaussian mean; should ==tv
  double mu = tv; 
  //need to be careful with sd, can get -ve values otherwise
  //this needs to be small so that last shell injection time
  //is ~tw
  double sd = 0.001; //standard deviation
  
  Numerical nums;
  int idum= (-1*time(NULL));
   
    
  double random1 = 0.0, random2 = 0.0, random3 = 0.0;
  double dt = 0.0; //injection interval
  double m = 0.0;
  //subtracting 1. for the puposes of generating a random 
  //number betwenn 0 and gmax. 1 is added after the obtaining 
  //the random number (e.g. for g between 1 and 1.9, generate
  //random number between 0 and 0.9 and then add one to it).
  double gMin = extGamMin;
  double g_max = (extGamMax - gMin);
    
  double gam = 0.0;
  //double shWidth = extShWidth; //Shell width
  double shWidth = 0.0;
  double gap=extShWidth;
  double beta=0.0;
  
  if(gap >= 1.0)
        {
          std::cerr<<"****Warning! Decrease width factor to less than 1.0 *****"<<"\n"
                   <<"****Will cause pileup at the source and errors*****"<<"\n"
                   <<"Exiting the program"
                   <<std::endl;
          exit(0);
          
        }
  
  std::vector<double>::iterator injTime, ma, ga, le;
  std::string place1("populateWithShells; useShellFile; Container.cc");
    
  
  try
    {
      
      std::string filename = extShellFile;
      std::ofstream shFile(filename.c_str());
      
      if (shFile.is_open())
        {
          shFile<<std::setiosflags(std::ios::scientific)<<std::setprecision(9);
          std::cout<<"Over-writing the file \""<<filename
                   <<"\" with randomly generated values"<<std::endl;              
          for(int i=0; i<=n; ++i)
            { 
              //standard rand() is just aweful!!
                  
              //avoid gamma = 0.0, adding 1 makes sure it is above 1 i.e. moving
                  
              random1=nums.random(idum);
              
              random2=nums.random(idum);
              
              random3=nums.random(idum);
              
              //gMu = (gamMax - gamMin);
              gam = gMin+(random3*g_max);
              //gam = (gaussian(gMu, sdGam, random3, random4));
              beta = sqrt(1. - (1. / pow(gam, 2)));
              //can use round(double) function if int required
              //need the time interval dist to be gaussian NOT inj t
              dt = (gaussian(mu, sd, random1, random2));     
                            
              m = (lum * dt) / (gam * sqr(physcon.c));
              
              shWidth = beta*physcon.c*dt*gap;
                            
              ti += dt;
              
              shFile<<ti<<"\t"<<m<<"\t"<<gam<<"\t"<<shWidth<<std::endl;
              
              
            }
        }
      else 
        {
          std::cout<<"Error! Unable to open shell file to write"<<std::endl;
        }
      
      
      
      
    }
  catch(std::bad_alloc) 
    {
      std::cerr<<"Memory exhauseted in populateWithShells;!useShellFile;Container.cc\n";
    }
  

  
}
//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
void Container::updateShellParams(double &t)
{
  bool inVacuum = extInVacuum;
  bool slowEnergization = extSlowEnerg;
  bool injWithIntEnergy = extInjIntEne;
  for (iterator sh = this->begin();
       sh != this->end(); sh++)
    {
      
      //(*sh)->radiativeLosses(t);
         
      if(!inVacuum)
        { 
	  (*sh)->adiabaticLosses(t);
	}
      
      (*sh)->setDimensions(t);
      
      (*sh)->setExpansionBeta();
   
      if(slowEnergization)
        {
          (*sh)->slowEnergization(t);
        }
            
      if(inVacuum)
        {
          (*sh)->inVacuumExpansion();
        }
     
      if(injWithIntEnergy)
        {
          (*sh)->shockZone(t);
        }
            
    }
}

//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
void Container::shellData(std::vector<double> &inj, std::vector<double> &ma,
                          std::vector<double> &ga, std::vector<double> &len)
{
  try
    {
      
      double time, mass, gamma, length;
      
      std::string filename = extShellFile;
      std::ifstream datafile(filename.c_str());
      std::cout<<"Reading shells parameters file \""<<filename<<"\"...";
      if (!datafile.good())
        throw std::runtime_error("Data file failed to load");
      
      while(!datafile.eof())
        {
          std::string line;
          getline(datafile, line);
          std::stringstream stream;
          
          if(line != "")
            {
              stream << line;
              stream >> time >> mass >> gamma >> length;
              //std::cout<<time<<"\t"<<mass<<"\t"<<gamma<<"\t"<<length<<std::endl;
              
              inj.push_back(time);
              ma.push_back(mass);
              ga.push_back(gamma);
              len.push_back(length);
            }
          
        }
      
      
      datafile.close();
       std::cout<<" done"<<std::endl;
    }
  catch(std::bad_alloc)
    {
      std::cerr<<"Memory exhauseted in shellData(arguments); Container.cc\n";
    }
}

//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
double Container::gaussian(const double &mean, const double &std,
                             double r1 , double r2)
{ 
    
  const double pi=3.1415927;
  //Box-Muller transformation method
  return (((sqrt(-2*log(r1)) * cos(2*pi*r2)) * std) + mean);

} 

//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
void Container::coutContainerElements()
{
for (const_iterator j = this->begin();
       j != this->end(); ++j)
    {  
      std::cout<<"all the elements"<<"\t"<<(*j)->getTimeOfInjection()<<"\t"
               <<(*j)->getMergeCount()<<"\t"
               <<(*j)->getOuterRadius()<<"\t"
               <<(*j)->getInnerRadius()<<"\t"
               <<(*j)->getShellWidth()<<"\t"
               <<(*j)->getInternalEnergy()<<"\n";
      
        //<<(*j)->getMergeCount()<<"\n";
      
      
               
               
      
    }
}
//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
void Container::majorEjection(double &tiM, double &tNow)
{
  try
    {
      std::cout<<"Major ejection event at (s): "<<tiM<<std::endl;
      Shells *lastShell;
      double lum = extJetLum;
      double gamM = 3.;
      double shWidthM = 1.e10;
      double mM = (lum * (tiM-tNow)) / (gamM * sqr(physcon.c));
      double inEneM = 0.9 * gamM * mM * (physcon.c * physcon.c);
              
      this->push_back(new Shells(tiM, mM, gamM, inEneM, shWidthM));
      lastShell = this->back();
      lastShell->setIdFalse();
      lastShell->setDimensions(tiM);
      // lastShell->setLowerAndUpperX(innerRadius, 
      //                                 outerRadius);
      //   lastShell->setInitialMergedVA();
      lastShell->powerlawNorm();
      lastShell->setExpansionBeta();
      lastShell->initialBEnergyDensity();
    }
  catch(std::bad_alloc)
    {
      std::cerr<<"Memory exhauseted in getData; Container.cc\n";
    }
}
//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
void Container::overlapPrevention(double &tNow)
{
  //the following picks up on any overlapping shells and artificially
  //fixes the problem. This can be avoided as long as 
  //(shell_length / injection gap) <= 1e6
  
  iterator inner, outer;
  double innerRad, outerRad;
  double setRad;
  Shells *innerShell;
  Shells *outerShell;
  
  outer = this->begin();
  inner = this->begin();
  ++inner;
  
  while(inner != this->end())
    {
      outerShell = *(outer);
      innerShell = *(inner);
      
      outerRad = outerShell->getInnerRadius();
      innerRad = innerShell->getOuterRadius();
      double t1 = innerShell->getTimeOfInjection();
      double t2 = outerShell->getTimeOfInjection();
            
      if(innerRad >= outerRad)
        {
          //std::cout<<"Overlap occured at ~ "<<"\t"<<outerRad<<"\t"<<innerRad<<"\t"<<t1<<"\t"<<t2<<std::endl;
          std::cerr<<"......................................................................."<<"\n";
          std::cerr<<"Overlap/pile up occured at ~ "<<"\t"<<outerRad<<"\t"<<"m"<<"\n";
          std::cerr<<"Time at pile up (s):"<<"\t"<<tNow<<"\n";
          std::cerr<<"Fixed for shells (inj t):"<<"\t"<<t1<<"\t"<<t2
                   <<"\n"<<"Reduce shell_width_factor (< 1.0); increase dt_inj"
                   <<"\n"<<"Increase BLF_min"<<std::endl;
          setRad = outerRad - 10.;
          innerShell->setOuterRadius(setRad);
        }
      ++outer;
      ++inner;
      
    }
  //std::cout<<"end"<<std::endl;
  
 
     
}
//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
bool Container::fileExists(const std::string& fileName)
{
  std::fstream fin;
  fin.open(fileName.c_str(),std::ios::in);
  if( fin.is_open() )
    {
      fin.close();
      return true;
    }
  fin.close();
  return false;
}
//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
void Container::injectShell(double &time, double &mass, double &gamma, 
                            double &length)
{
  double inEne = 0.0;
  //std::auto_ptr<Shells> p(new Shells(time, mass, gamma, inEne, length));
  this->push_back(new Shells(time, mass, gamma, inEne, length));
  //this->push_back(p);
  Shells *newShell = this->back();
  
  newShell->shellInitialization();
  
                  
}
//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
void Container::dumpActiveShells()
{
  try
    {
      
      std::ofstream activeFile("activeShells.dat");
      std::cout<<"Writing all the active shell feature to \"activeShells.dat\"....";
      if (!activeFile.good())
        throw std::runtime_error("activeShell.dat file failed to open");
      
      Container::const_iterator shells;
      double ti, mass, gamma, ienergy, width, betaExpansion, x_l, x_u;
      double r_l, r_u, vol, area, egammaMin, egammaMax, pNorm, BEDens; 
      double fracEInt, intEavail;
      double magPressure, shellCentre, tPrev, shellThermEner, tEnergize;
      int mergeCount;
      bool shell, shocked;
      std::vector<double> iNuGrid;
      std::vector<double> tau;
      std::vector<double>::iterator iNuIt, tauIt;
      
      activeFile<<std::setiosflags(std::ios::scientific)<<std::setprecision(8);
      
      for (shells = this->begin(); shells != this->end(); ++shells)
        {
          ti = (*shells)->getTimeOfInjection();
          mass = (*shells)->getShellMass();
          gamma = (*shells)->getShellGamma();
          ienergy = (*shells)->getInternalEnergy();
          width = (*shells)->getShellWidth();
          betaExpansion = (*shells)->getExpansionBeta();
          x_l = (*shells)->getInnerRadius();
          x_u = (*shells)->getOuterRadius();
          r_l = (*shells)->getShellInRad();
          r_u = (*shells)->getShellOuRad();
          vol = (*shells)->getShellVolume();
          area = (*shells)->getShellArea();
          egammaMin = (*shells)->getGammaMin();
          egammaMax = (*shells)->getGammaMax();
          pNorm = (*shells)->getPLawNorm();
          BEDens = (*shells)->getBEneDens();
          fracEInt = (*shells)->getIntEfrac();
          intEavail = (*shells)->getIntEavail();
          magPressure = (*shells)->getMagPressure();
          shellCentre = (*shells)->getLocation();
          tPrev = (*shells)->getPrevTime();
          shellThermEner = (*shells)->getThermalEnergy();
          tEnergize = (*shells)->getTEnergize();
          mergeCount = (*shells)->getMergeCount();
          shell = (*shells)->getShellId();
          shocked = (*shells)->shellShocked();
          iNuGrid = (*shells)->iNuVals();
          tau = (*shells)->tauVals();
          
          activeFile<<"Shell"<<"\t";
          activeFile<<ti<<"\t"<<mass<<"\t"<<gamma<<"\t"<<ienergy<<"\t"
                      <<width<<"\t"<<betaExpansion<<"\t"<<x_l<<"\t"<<x_u
                      <<"\t"<<r_l<<"\t"<<r_u<<"\t"<<vol<<"\t"<<area<<"\t"
                      <<egammaMin<<"\t"<<egammaMax<<"\t"<<pNorm<<"\t"
                      <<BEDens<<"\t"<<fracEInt<<"\t"<<intEavail<<"\t"
                      <<magPressure<<"\t"<<shellCentre<<"\t"<<tPrev<<"\t"
                      <<shellThermEner<<"\t"<<tEnergize<<"\t"<<mergeCount<<"\t"
                      <<shell<<"\t"<<shocked<<"\n";
          activeFile<<"Nu"<<"\t";
          for(iNuIt = iNuGrid.begin(); iNuIt != iNuGrid.end(); ++iNuIt)
            {
              activeFile<<*iNuIt<<"\t";
            }
          activeFile<<"\n"<<"Tau"<<"\t";
          for(tauIt = tau.begin(); tauIt != tau.end(); ++tauIt)
            {
              activeFile<<*tauIt<<"\t";
            }
          
          activeFile<<"\n"<<"Inject"<<std::endl;
          
        }
      
      
      activeFile.close();
       std::cout<<" done"<<std::endl;
       std::cout<<"Number of shells written to file: "<<this->size()<<std::endl;
       
    }
  catch(std::bad_alloc)
    {
      std::cerr<<"Memory exhauseted in dumpActiveShell(); Container.cc\n";
    }
}
//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
void Container::readActiveShells()
{
  try
    {
      
      std::ifstream activeFile("activeShells.dat");
      std::cout<<"Reading the previous active shell feature from \"activeShells.dat\"....";
      if (!activeFile.good())
        throw std::runtime_error("activeShell.dat file failed to open");
      
      Container::const_iterator shells;
      double ti, mass, gamma, ienergy, width, betaExpansion, x_l, x_u;
      double r_l, r_u, vol, area, egammaMin, egammaMax, pNorm, BEDens; 
      double fracEInt, intEavail;
      double magPressure, shellCentre, tPrev, shellThermEner, tEnergize;
      int mergeCount;
      bool shell, shocked;
      std::vector<double> iNuGrid;
      std::vector<double> tau;
      std::vector<double>::iterator iNuIt, tauIt;
      std::string temp;
      double temp2, temp3;
      
      while(!activeFile.eof())
        {
          std::string line;
          getline(activeFile, line);
          std::stringstream stream;
          
          if(line != "")
            {
              stream << line;
              stream >> temp;
              if(temp == "Shell")
                {
                  stream >> ti >> mass >> gamma >> ienergy >> width >> betaExpansion >> x_l 
                     >> x_u >> r_l >> r_u >> vol >> area >> egammaMin >> egammaMax >> pNorm
                     >> BEDens >> fracEInt >> intEavail >> magPressure >> shellCentre 
                     >> tPrev >> shellThermEner >> tEnergize >> mergeCount >> shell 
                     >> shocked; 
                }
              if(temp == "Nu")
                {
                  while(!stream.eof())
                    {
                      stream >> temp2;
                      iNuGrid.push_back(temp2);
                    }
                }
              if(temp == "Tau")
                {
                 while(!stream.eof())
                    {
                      stream >> temp3;
                      tau.push_back(temp3);
                    }
                } 
              if(temp == "Inject")
                {
                  
              this->push_back(new Shells(ti, mass, gamma, ienergy, width));
              Shells *newShell = this->back();
              newShell->setExpansionBeta(betaExpansion);
              newShell->setInnerRadius(x_l);
              newShell->setOuterRadius(x_u);
              newShell->setShellInRad(r_l);
              newShell->setShellOuRad(r_u);
              newShell->setShellVolume(vol);
              newShell->setShellArea(area);
              newShell->setGammaMin(egammaMin);
              newShell->setGammaMax(egammaMax);
              newShell->setPLawNorm(pNorm);
              newShell->setBEneDens(BEDens);
              newShell->setIntEfrac(fracEInt);
              newShell->setIntEavail(intEavail);
              newShell->setMagneticPressure(BEDens);
              newShell->setLocation(shellCentre);
              newShell->setPrevTime(tPrev);
              newShell->setThermalEnergy(shellThermEner);
              newShell->setTEnergize(tEnergize);
              newShell->setMergeCount(mergeCount);
              newShell->setShellId(shell);
              newShell->setShellShocked(shocked);
              newShell->copyINu(iNuGrid);
              newShell->copyTau(tau);
              //std::cout<<"test "<<ti<<std::endl;
                }
              
            }
          
        }
      
      activeFile.close();
       std::cout<<" done"<<std::endl;
     

    }
  catch(std::bad_alloc)
    {
      std::cerr<<"Memory exhauseted in readActiveShell(); Container.cc\n";
    }
}
//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
// void Container::getData()
// {
//   try
//     {
      
//   timeAxis = new std::vector<double>;
//   fluxAxis = new std::vector<double>;
//   double timeData, fluxData, errorData;
  
//   std::string filename = extDatFile;
//   std::string line;
//   std::ifstream datafile(filename.c_str());
//   if (!datafile.good())
//     throw std::runtime_error("Data file failed to load");
//   std::cout<<"Reading the lightcurve..."<<std::endl;
  
//   while(!datafile.eof())
//     {
//       datafile >> timeData >> fluxData >> errorData;
//       timeAxis->push_back(timeData);
//       fluxAxis->push_back(fluxData);
            
//     }
  
//   datafile.close();

//   std::string shellFile = extShellFile;
//   //const std::string shFileCheck = extShellFile;
//   std::ofstream shellData(shellFile.c_str(), std::ofstream::app);
//   std::vector<double>::iterator t, f;
//   double dt=0.0, g_mean = 0.0, tot=0.0, lum=0.0, random1=0.0, random2=0.0;
//   double ti=1.0, gam=0.0, dl=0.0;
//   double mass=1.e12;
//   double sd = 0.01;
  
//   Numerical numer;
//   int idum= (-1*time(NULL));
//   //bool shFileExist = fileExists(shFileCheck);
  
//   //if(!shFileExist)
//   //  {
//   if (shellData.is_open())
//     {
//       std::cout<<"Writing shells parameters file..."<<std::endl;
          
//       for(t=timeAxis->begin(), f=fluxAxis->begin();
//           t != timeAxis->end(); ++t, ++f)
//         {
//           lum = *f*1e30;
//           if(lum <= 1.e29)
//             {
//               g_mean=1.6;
//             }
//           else if(lum >= 1.e29 && lum <= 1.e30)
//             {
//               g_mean=1.7;
//             }
//           else if(lum >= 1.e30 && lum <= 1.2e30)
//             {
//               g_mean=1.8;
//             }
//           else if(lum >= 1.2e30 && lum <= 1.4e30)
//             {
//               g_mean=1.9;
//             }
//           else if(lum >= 1.4e30 && lum <= 1.6e30)
//             {
//               g_mean=2.0;
//             }
          
          
//           dt=(g_mean*mass*sqr(physcon.c)) / lum;
//           tot = 0.0;
//           //0.5 added as the time starts at 0.5
//           //and by starting at 1. and then carrying, avoids the
//           //overlap you would get otherwise.
//           //if the t starts at 1. then remove this 0.5 factor
//           ti=(*t)+0.5; 
//           //std::cout<<*t<<"\t"<<ti<<std::endl;
          
//           while(tot <= 0.9)
//             {
//               random1=numer.random(idum);
//               random2=numer.random(idum);
//               gam=(gaussian(g_mean, sd, random1, random2));
//               //gam = 1.+(random1*(gmax-1.));
              
//               if (gam == 1.){gam+=0.0001;}
              
//               //std::cout<<gam<<"\t"<<random1<<std::endl;
//               ti += dt;
//               tot += dt;
//               dl=(1.e5)*dt;
//               //std::cout<<ti<<"\t"<<dt<<std::endl;
//               shellData<<ti<<"\t"<<mass<<"\t"<<gam<<"\t"<<dl<<"\n";
//             }
//           ti=0.0;
//           dt=0.0;
          
          
//         }
      
    
//     }
//   else 
//     {
//       std::cerr<<"Error! Unable to open write shell data file"<<std::endl;
//     }
  
//   shellData.close();
//   //  }
//   // else 
//   //   {
//   //    break;
//   //   }
//     }
//   catch(std::bad_alloc)
//     {
//        std::cerr<<"Memory exhauseted in getData; Container.cc\n";
//     }
  
// }
// //.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
// void Container::shellData()
// {
//   try
//     {
      
//       shITime = new std::vector<double>;
//       shM = new std::vector<double>;
//       shGam = new std::vector<double>;
//       shLength = new std::vector<double>;
//       double time, mass, gamma, length;
      
//       std::string filename = extShellFile;
//       std::ifstream datafile(filename.c_str());
//       std::cout<<"Reading shells parameters file \""<<filename<<"\"...";
//       if (!datafile.good())
//         throw std::runtime_error("Data file failed to load");
      
//       while(!datafile.eof())
//         {
//           std::string line;
//           getline(datafile, line);
//           std::stringstream stream;
          
//           if(line != "")
//             {
//               stream << line;
//               stream >> time >> mass >> gamma >> length;
              
//               shITime->push_back(time);
//               shM->push_back(mass);
//               shGam->push_back(gamma);
//               shLength->push_back(length);
//             }
          
//         }
      
//       datafile.close();
//       std::cout<<" done"<<std::endl;
      
//     }
//   catch(std::bad_alloc)
//     {
//       std::cerr<<"Memory exhauseted in shellData; Container.cc\n";
//     }
// }
//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
