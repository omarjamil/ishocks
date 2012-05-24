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

#include <iostream>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <algorithm>
#include <string>

#include "evolution.hh"
#include "funcobj.hh"
#include "physcon.hh"
#include "paraccess.hh"
#include "radiation.hh"

//used for multi runs and search through parameter space
//also produces a file with average spectral index and spectral index
//in the lightcurve file. Good idea to only have individual frequencies when
//using GRID_SEARCH
//#define GRID_SEARCH

//****TEMPORARY??****
//reads the file splitLC.par
#define SPLIT_LC
//****??TEMPORARY??****

//if TAU_WRITE defined then the finalTStau.dat file is written. This file
//has the optical depths of all the shells at final time step.
#define TAU_WRITE

//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....



// evolution class constructor
// all and active initialized to the Shell deques passed
Evolution::Evolution(Container *active_)
  :activeShells(active_)
{}

//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....

// evolution class destructor
Evolution::~Evolution()
{
  delete activeShells;
  delete nu;
  delete fluxNu1;
  delete fluxNu2;
}

//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....

void Evolution::evolve(int &duration)
{
  try
    {

      Container::const_iterator all, next;
      Container::iterator locator, active, inLocate, ouLocate;

      double col_t;
      double dt_collision, dt_next, t_now, t_prev;
      double t_inner = 0.0, t_outer = 0.0;
      double t_last;
      double t_next;
      double t_incr;
      const std::string exitFile = "exitNow";
      bool exitF;
      bool readPrev = extLoadPrevSim;
      bool dumpAtEnd = extDumpActive;

      //the range of freq for getting the synchrtron spectrum
      nu = new std::vector<double>;
      egamma = new std::vector<double>;
      
      //populating the gamma container (duplicate of the one in Shells)
      gammaRange(egamma);
      
      //populating the freq container
      freqRange(nu);

      //Compton Table
      compTable(nu, egamma);

      //For getting the average flux
      fluxNu1 = new std::vector<double>;
      fluxNu2 = new std::vector<double>;
      //the following vector contain shell inj time, mass, gamma, length
      //that are used to inject new ones.
      std::vector<double> shellI;
      std::vector<double> shellM;
      std::vector<double> shellG;
      std::vector<double> shellL;
      std::vector<double>::const_iterator shellIT, shellITNext;
      std::vector<double>::const_iterator shellMT;
      std::vector<double>::const_iterator shellGT;
      std::vector<double>::const_iterator shellLT;
      int tCountPrev = 0;

      //load the previous simulation
      if(readPrev)
        {
          activeShells->readActiveShells();
          tCountPrev = activeShells->size();
        }
      // does not do anything to activeShells
      // just a way to call shellData routine here
      //activeShells->coutContainerElements();
      activeShells->shellData(shellI, shellM, shellG, shellL);

      //Radiation object to obtain the spectrum
      Radiation spectrum;

      //contains the merger shells info
      //returned by minCollisionTime
      std::vector<Mergers*> mergerEvent;
      //the last shell
      t_last = shellI.back();

      double iS, mS, gS, lS;
      //the very first shell; will not collide with anything
      shellIT = shellI.begin();
      shellMT = shellM.begin();
      shellGT = shellG.begin();
      shellLT = shellL.begin();
      iS = *shellIT; mS = *shellMT; gS = *shellGT; lS = *shellLT;
      activeShells->injectShell(iS, mS, gS, lS);

      //the next shell
      shellITNext = shellI.begin();

      //t_now = (*all)->getTimeOfInjection();
      t_now = *shellIT;


      //One shell ahead of "all"
      ++shellITNext;
      t_next = *shellITNext;

      int tCount = 1;
      int tCount2 = 1;
      int mCount = 0;

      std::cout<<"Running the simulation, jet on (shells + mergers)..."<<std::endl;

      while(*shellITNext != shellI.back())
        {
          //the shell to be injected into activeshells
          ++shellIT; ++shellMT; ++shellGT; ++shellLT;
          iS = *shellIT; mS = *shellMT; gS = *shellGT; lS = *shellLT;
          t_now = *shellIT;

          activeShells->injectShell(iS, mS, gS, lS);
          activeShells->updateShellParams(t_now);

          //to check for overlaps
          activeShells->overlapPrevention(t_now);

          spectrum.synchSpectrum(nu, activeShells);
          //std::cout<<t_now<<"\t"<<"injection"<<std::endl;

#ifndef GRID_SEARCH
          results(activeShells, nu, t_now, tCount);
#endif

#ifdef GRID_SEARCH
          results(activeShells, nu, t_now, fluxNu1, fluxNu2);
#endif
          //time step count
          tCount += 1;
          //tCount2 += 1;

          //std::cout<<"\r"<<"Time step.................................: "<<tCount2;

          //the next shell to be injected
          ++shellITNext;
          t_next = *shellITNext;
          dt_next = t_next - t_now;

          //the following loop does mergers as long as they are before the
          //time for next injection
          do
            {
              mergerEvent.clear();
              dt_collision = minCollisionTime(activeShells, mergerEvent);

              if (dt_collision > 0.0 && dt_collision <= dt_next)
                {
                  if(!mergerEvent.empty())
                    {

                      t_now += dt_collision;
                      dt_next -= dt_collision;

                      merger = new Mergers(*mergerEvent[0]);
                      //deleting the mergers objects and pointers.
                      std::for_each(mergerEvent.begin(), mergerEvent.end(), delete_ptr<Mergers>());
                      mergerEvent.clear();

                      activeShells->updateShellParams(t_now);
                      merger->doMerger(activeShells, t_now);
                      delete merger;
                      spectrum.synchSpectrum(nu, activeShells);
                      //std::cout<<t_now<<"\t"<<"merger"<<std::endl;
#ifndef GRID_SEARCH
                      results(activeShells, nu, t_now, tCount);
#endif

#ifdef GRID_SEARCH
                      results(activeShells, nu, t_now, fluxNu1, fluxNu2);
#endif
                      activeShells->overlapPrevention(t_now);
                      mCount += 1;
                      //tCount2 += 1;
                      //std::cout<<"\r"<<"Time step.................................: "<<tCount2;
                      //to catch any collision for the newly merged shell
                      dt_collision = minCollisionTime(activeShells, mergerEvent);
                    } else
                    {
                      std::cerr<<"Error! No Merger calculated"<<"\n";
                    }
                }
              else
                {
                  break;
                }
              //dt_collision = minCollisionTime(activeShells, mergerEvent);
            }
          while (dt_collision <= dt_next);
          //std::cout<<t_now<<std::endl;

          //this should inject the last shell
          if (*shellITNext == t_last)
            {
              iS = shellI.back();
              mS = shellM.back();
              gS = shellG.back();
              lS = shellL.back();
              activeShells->injectShell(iS, mS, gS, lS);
              tCount += 1;
              //std::cout<<gS<<std::endl;

              break; //to break the while loop at the end of allShells
            }

          //this checks the if "exitNow" file exists. If it is created then
          //the programme exits and writes the final time step before doing so.
          exitF = fileExists(exitFile);
          if(exitF)
            {
              std::cout<<"exitNow file detected"<<"\n"
                       <<"writing final time step files"<<std::endl;


              bool finalTS = extFinalTStep;
              if(finalTS)
                {

                  finalTStep(activeShells, nu, t_now);
#ifdef TAU_WRITE
                  finalTStepTau(activeShells, nu, t_now);
#endif
                }
              std::cout<<"calling exit(0) in Evolution::evolve"<<std::endl;
              activeShells->dumpActiveShells();
              exit(0);

            }


        }
      //to write the last time step to the "lastTStep.dat"
      bool finalTS = extFinalTStep;
      if(finalTS)
        {

          finalTStep(activeShells, nu, t_now);

#ifdef TAU_WRITE
          finalTStepTau(activeShells, nu, t_now);
#endif
        }

      //To continue evolution after jet has been turned off...
      //continueEvol(t_now, activeShells);

#ifdef GRID_SEARCH
      avgLC(fluxNu1, fluxNu2);
#endif
      std::cout<<"Number of shells loaded from previous sim: "<<tCountPrev<<"\n";
      std::cout<<"Number of shells injected: "<<tCount<<"\n";
      std::cout<<"Number of mergers: "<<mCount<<"\n";
      std::cout<<"Number of time steps: "<<tCount + mCount<<"\n";
      std::cout<<"Time since the start of the simulation(s) "<<t_now<<"\n";
      if(dumpAtEnd)
        {
          activeShells->dumpActiveShells();
        }
    }
  catch(std::bad_alloc)
    {
      std::cerr<<"Memory exhauseted in evolve;Evolution.cc\n";
      activeShells->dumpActiveShells();
    }

}


//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....

void Evolution::evolveIncrSampling(int &duration)
{
  try
    {

      std::cout<<"Increase Time sampling, making it \"more\" discrete steps"<<std::endl;
      std::cout<<"Time step resolution (~s): "<<extStepResolution<<std::endl;
      std::cout<<"Total run time (s): "<<extTotalDuration<<std::endl;

      Container::const_iterator all, next;
      Container::iterator locator, active, inLocate, ouLocate;

      const std::string exitFile = "exitNow";
      bool exitF;
      bool readPrev = extLoadPrevSim;
      bool dumpAtEnd = extDumpActive;

      double col_t, dt_collision, dt_next, t_now = 0.0;
      double t_last, t_next, dtIncr;
      double t_max = extTotalDuration;
      double t_incr = extStepResolution;
      int tCount = 1, mCount = 0, shCount = 0; //timestep and merger count
      int tCountPrev = 0;
      nu = new std::vector<double>; //container for the frequencies
      freqRange(nu); //populating the freq container
      //For getting the average flux
      fluxNu1 = new std::vector<double>;
      fluxNu2 = new std::vector<double>;
      //the following vector contain shell inj time, mass, gamma, length
      //that are used to inject new ones.
      std::vector<double> shellI;
      std::vector<double> shellM;
      std::vector<double> shellG;
      std::vector<double> shellL;
      std::vector<double>::const_iterator shellIT, shellITNext;
      std::vector<double>::const_iterator shellMT;
      std::vector<double>::const_iterator shellGT;
      std::vector<double>::const_iterator shellLT;

      //load the previous simulation
      if(readPrev)
        {
          activeShells->readActiveShells();
          tCountPrev = activeShells->size();
        }

      // does not do anything to activeShells
      // just a way to call shellData routine here
      activeShells->shellData(shellI, shellM, shellG, shellL);

      //Ghost shell necessary for injecting the final shell
      //the ghost shell itself is not injected.
      double tGhost = 1.e-6; double mGhost = 1.; double gGhost = 1.;
      double lGhost = 1.;
      shellI.push_back(tGhost);
      shellM.push_back(mGhost);
      shellG.push_back(gGhost);
      shellL.push_back(lGhost);

      Radiation spectrum; //Radiation object to obtain the spectrum
      std::vector<Mergers*> mergerEvent; //contains the merger shells info
      //returned by minCollisionTime

      //the final shell time
      t_last = shellI.back();


      //the very first shell;
      //will not collide with anything

      double iS, mS, gS, lS;
      shellIT = shellI.begin();
      shellMT = shellM.begin();
      shellGT = shellG.begin();
      shellLT = shellL.begin();
      iS = *shellIT; mS = *shellMT; gS = *shellGT; lS = *shellLT;

      activeShells->injectShell(iS, mS, gS, lS);

      shCount += 1;
      t_now = *shellIT;

      //the next shell
      shellITNext = shellI.begin();

      if(*(shellIT) != shellI.back())
        {
          //One shell ahead of "all"
          ++shellITNext;
          t_next = *shellITNext;
          dt_next = t_next - t_now;

        }
      std::cout<<"Running the simulation, jet on (shells + mergers)..."<<std::endl;

      while(t_now < t_max)
        {


          dtIncr = t_incr;

          if(*(shellITNext) == t_last)
            {
              dt_next = t_max; //So that the conditions are met in order to
              //carry on collisions after the injections are complete.
            }
          else
            {
              dt_next = *shellITNext - t_now;
            }


          dt_collision = minCollisionTime(activeShells, mergerEvent);
          //std::cout<<dtIncr<<"\t"<<dt_next<<"\t"<<dt_collision<<std::endl;

          //Loop for injecting shells
          while(dt_next <= dtIncr && dt_next < dt_collision)
            {
              if(*(shellITNext) == t_last)
                {
                  break; //to break the while loop at the end of allShells
                }

              ++shellIT; ++shellMT; ++shellGT; ++shellLT;
              iS = *shellIT; mS = *shellMT; gS = *shellGT; lS = *shellLT;
              t_now = *shellIT;

              activeShells->injectShell(iS, mS, gS, lS);

              shCount += 1;
              activeShells->updateShellParams(t_now);
              activeShells->overlapPrevention(t_now); //to check for overlaps
              spectrum.synchSpectrum(nu, activeShells);//the spectrum from the shells.
#ifndef GRID_SEARCH
              results(activeShells, nu, t_now, tCount);
#endif

#ifdef GRID_SEARCH
              results(activeShells, nu, t_now, fluxNu1, fluxNu2);
#endif
              tCount += 1; //time step count

              if(*(shellIT) != shellI.back())
                {
                  //One shell ahead of "all"
                  ++shellITNext;
                  t_next = *shellITNext;
                }

              dtIncr -= dt_next;
              dt_collision = minCollisionTime(activeShells, mergerEvent);

              t_next = *shellITNext;
              dt_next = t_next - t_now;

            }
          //Loop for shell collisions
          while(dt_collision < dt_next && dt_collision < dtIncr)
            {
              mergerEvent.clear();
              dt_collision = minCollisionTime(activeShells, mergerEvent);
              if (dt_collision <= dt_next && dt_collision <dtIncr)
                {
                  if(!mergerEvent.empty())
                    {
                      t_now += dt_collision;

                      dt_next -= dt_collision;
                      dtIncr -=dt_collision;
                      merger = new Mergers(*mergerEvent[0]);
                      //deleting the mergers objects and pointers.
                      std::for_each(mergerEvent.begin(), mergerEvent.end(), delete_ptr<Mergers>());
                      mergerEvent.clear();
                      activeShells->updateShellParams(t_now);
                      merger->doMerger(activeShells, t_now);
                      delete merger;

                      //to catch any collision for the newly merged shell
                      dt_collision = minCollisionTime(activeShells, mergerEvent);
                      spectrum.synchSpectrum(nu, activeShells);
#ifndef GRID_SEARCH
                      results(activeShells, nu, t_now, tCount);
#endif

#ifdef GRID_SEARCH
                      results(activeShells, nu, t_now, fluxNu1, fluxNu2);
#endif
                      activeShells->overlapPrevention(t_now);
                      mCount += 1;
                      tCount += 1;

                      //std::cout<<"\r"<<"Time step.................................: "<<tCount2;
                    } else
                    {
                      std::cerr<<"Error! No Merger calculated"<<"\n";
                    }
                }
              else
                {
                  break;
                }
            }

          if(dtIncr > 0.0 && dtIncr < dt_next && dtIncr <dt_collision)
            {
              t_now += dtIncr;
              activeShells->updateShellParams(t_now);
              activeShells->overlapPrevention(t_now);
              spectrum.synchSpectrum(nu, activeShells);

#ifndef GRID_SEARCH
              results(activeShells, nu, t_now, tCount);
#endif

#ifdef GRID_SEARCH
              results(activeShells, nu, t_now, fluxNu1, fluxNu2);
#endif
              tCount += 1;
            }

          //this checks the if "exitNow" file exists. If it is created then
          //the programme exits and writes the final time step before doing so.
          exitF = fileExists(exitFile);
          if(exitF)
            {
              std::cout<<"exitNow file detected"<<"\n"
                       <<"writing final time step files"<<std::endl;


              bool finalTS = extFinalTStep;
              if(finalTS)
                {

                  finalTStep(activeShells, nu, t_now);
#ifdef TAU_WRITE
                  finalTStepTau(activeShells, nu, t_now);
#endif
                }
              std::cout<<"calling exit(0) in Evolution::evolutionIncrSampling"<<std::endl;
              activeShells->dumpActiveShells();
              exit(0);

            }

        }

      //to write the last time step to the "lastTStep.dat"
      bool finalTS = extFinalTStep;
      if(finalTS)
        {
         finalTStep(activeShells, nu, t_now);

#ifdef TAU_WRITE
          finalTStepTau(activeShells, nu, t_now);
#endif
        }

#ifdef GRID_SEARCH
      avgLC(fluxNu1, fluxNu2);
#endif
      std::cout<<"Number of shells loaded from previous sim: "<<tCountPrev<<"\n";
      std::cout<<"Number of shells injected: "<<shCount<<"\n";
      std::cout<<"Number of mergers: "<<mCount<<"\n";
      std::cout<<"Number of time steps: "<<tCount + mCount<<"\n";
      std::cout<<"Time since the start of the simulation(s) "<<t_now<<"\n";
      if(dumpAtEnd)
        {
          activeShells->dumpActiveShells();
        }
    }
  catch(std::bad_alloc)
    {
      std::cerr<<"Memory exhauseted in evolveIncrSampling;Evolution.cc\n";
      activeShells->dumpActiveShells();
    }

}

//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....

double Evolution::minCollisionTime(Container *shells,
                                   std::vector<Mergers*> &merg)
{
  try
    {

      std::vector<Mergers*> mergingShells;
      std::vector<double> times;
      std::vector<double>::iterator t;
      Container::iterator inner, outer;

      double collTime, minTime;
      int itposition;
      double shI, shO, collTPrev=1.e8;

      outer = shells->begin();
      inner = shells->begin();
      ++inner;

      while(inner != shells->end())
        {
          collTime = collisionTime(*inner, *outer);

          if(collTime < collTPrev && collTime > 0.0)
            {
              collTPrev = collTime;
              shI = (*inner)->getTimeOfInjection();
              shO = (*outer)->getTimeOfInjection();
            }

          ++outer;
          ++inner;

        }

      if (collTPrev < 1.e8)
        {

          minTime = collTPrev;

          inner = std::find_if(shells->begin(), shells->end(),
                           FindShellWithDouble(shI));
          outer = std::find_if(shells->begin(), shells->end(),
                               FindShellWithDouble(shO));

          merg.push_back(new Mergers(*inner, *outer));

          return minTime;

        }
      else
        {
          return minTime = 1.e12;

        }

    }
  catch(std::bad_alloc)
    {
      std::cerr<<"Memory exhauseted in minCollisionTime;Evolution.cc\n";
    }

}
//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
double Evolution::collisionTime(Shells *inner, Shells *outer)
{
  double t_now = inner->getTimeOfInjection();
  double r_inner = inner->getOuterRadius();
  double r_outer = outer->getInnerRadius();

  //double r_inner = inner->getLocation();
  //double r_outer = outer->getLocation();

  double w_outer = outer->getShellWidth();
  double w_inner = inner->getShellWidth();
  double b_outer = outer->getShellBeta();
  double b_inner = inner->getShellBeta();
  double t_inj_inner = inner->getTimeOfInjection();
  double t_inj_outer = outer->getTimeOfInjection();
  double inner_beta_exp = inner->getExpansionBeta();
  double outer_beta_exp = outer->getExpansionBeta();

  //from Spada paper
  //this gives the time gap for shell collision
  // double t_col = ((r_inner - r_outer) - 0.5*(w_outer+w_inner))/
//     (((b_outer*physcon.c) - (b_inner*physcon.c)) +
//      0.5*(inner_beta_exp*physcon.c + outer_beta_exp*physcon.c));


//from my calculations; this avoids overlaps
// and better adapted for our model.
  double t_col = (r_outer - r_inner)/
    ((inner_beta_exp*physcon.c + outer_beta_exp*physcon.c)+
     ((b_inner*physcon.c) - (b_outer*physcon.c)));

  //This gives thea absolute collision time
 //  double t_col = ((b_inner * t_inj_inner - b_outer * t_inj_outer)/
//                   (b_inner - b_outer)) +((outer_beta_exp - inner_beta_exp)
//                                          /(physcon.c));


//  std::cout<<t_col<<"\t"<<r_outer<<"\t"<<r_inner<<"\t"<<std::endl;


  // if(t_col<0.0)
//     {
//       std::cout<<"t_col"<<"\t"<<t_col<<std::endl;
//     }

//  double t_col = (r_inner - r_outer) /
//    (((b_outer*physcon.c) - (b_inner*physcon.c)) +
//      0.5*(inner_beta_exp*physcon.c + outer_beta_exp*physcon.c));

//  double t_col = (r_outer - r_inner) /
//    (((b_inner*physcon.c) - (b_outer*physcon.c)) +
//     0.5*(inner_beta_exp*physcon.c + outer_beta_exp*physcon.c));


  return t_col;

}

//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
double Evolution::dtNextInjection(Shells *now,
                                      Shells *nextShell)
{
  double t_now = now->getTimeOfInjection();
  double t_next = nextShell->getTimeOfInjection();
  double dt = t_next - t_now;

  return dt;


}
//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
//will not produce any information regarding spectral index
//this will have to be done manually (post simulation) using the
//lightcurve file.
void Evolution::results(Container *activeSh, std::vector<double> *nu,
                        double &tNow, int &tCount)
{
  bool shell;
  double location;
  double ti;
  double shWidth;
  double innerRadius;
  double outerRadius;
  double bEneDens;
  double intEner;
  double shGamma;
  double shMass;
  double intEnerTot=0.0, totalEner=0.0;

  bool useResFile = extUseResFile;
  int every = extEveryTStep;

  std::string filename = extResultsFile;
  std::string lcFilename = extLCFile;
  std::ofstream timeFile(lcFilename.c_str(), std::ofstream::app);

  std::vector<double> nuVals(nu->size(), 0.0);

  std::vector<double>::iterator nuIt1;
  std::vector<double>::const_iterator nuValIt, nuItConst;

  //****TEMPORARY??****
  bool sp = extSplitLC;
  if(sp)
    {
      splitLightcurve(activeSh, tNow);
    }//****??TEMPORARY**** //****TEMPORARY??****
  else
    {
      //****??TEMPORARY****
  if(timeFile.is_open())
    {
      timeFile<<std::setiosflags(std::ios::scientific)<<std::setprecision(8)<<tNow;

      //timeFile<<tNow;

      for(Container::const_iterator sh = activeSh->begin();
          sh != activeSh->end(); ++sh)
          {
            shell = (*sh)->getShellId();
            std::vector<double> iNu = (*sh)->iNuVals();

            if(!shell)
              {
                for(nuValIt = iNu.begin(), nuIt1 = nuVals.begin();
                    nuValIt != iNu.end(); ++nuIt1, ++nuValIt)
                  {
                    *nuIt1 += *nuValIt/1.e-26;
                  }
              }
          }

        for(nuItConst = nuVals.begin();
            nuItConst != nuVals.end(); ++nuItConst)
          {
            timeFile<<"\t"<<*nuItConst;
          }

        timeFile<<"\n";
        nuVals.clear();

      }
  else
    {
      std::cerr<<"Error! Unable to open lightcurve file"<<std::endl;
    }
  timeFile.close();
  //****TEMPORARY??****
    }
  //****??TEMPORARY****

if(useResFile)
  {
    if(tCount % every == 0)
      {
        std::ofstream resultsFile(filename.c_str(), std::ofstream::app);

        if (resultsFile.is_open())
              {
                resultsFile<<std::setiosflags(std::ios::scientific)<<std::setprecision(8);

                for(Container::const_iterator sh = activeSh->begin();
                    sh != activeSh->end(); ++sh)
                  {
                    shell = (*sh)->getShellId();
                    location = (*sh)->getLocation();
                    ti = (*sh)->getTimeOfInjection();
                    innerRadius = (*sh)->getInnerRadius();
                    outerRadius = (*sh)->getOuterRadius();
                    shWidth = outerRadius - innerRadius;//(*sh)->getShellWidth();
                    bEneDens = (*sh)->getBEneDens();
                    intEner = (*sh)->getInternalEnergy();
                    shGamma = (*sh)->getShellGamma();
                    shMass = (*sh)->getShellMass();

                    totalEner += shGamma * shMass * physcon.c * physcon.c;
                    intEnerTot += intEner;

                    //std::cout<<bEneDens<<std::endl;

                    std::vector<double> iNu = (*sh)->iNuVals();

                    if(location >= 0.0 && location <= 1.e20)
                      {
                        resultsFile<<"\n"<<tNow<<"\t"<<location;
                      if(!shell)
                        {
                          for(std::vector<double>::const_iterator it = iNu.begin(),
                                it3= nu->begin(); it != iNu.end(); ++it, ++it3)
                            {
                              resultsFile<<"\t"<<*it/1.e-26;
                            }
                          resultsFile<<"\t"<<shWidth;
                        }
                      else
                        {
                          for(std::vector<double>::const_iterator it2 = nu->begin();
                              it2 != nu->end(); ++it2)
                            {
                              resultsFile<<"\t"<<0.0;

                            }
                          resultsFile<<"\t"<<shWidth;
                        }
                      }
                  }

                totalEner=0.0;
                intEnerTot=0.0;

              resultsFile<<"\n";
              }
        else
          {
            std::cout<<"Error! Unable to open results file"<<std::endl;
          }
        resultsFile.close();

        //this writes a file that contains every time step with
        //with each shell location and optical depths. The optical depths
        //will correspond to the frequencies given in nuRange.dat
#ifdef TAU_WRITE
        tauWrite(activeShells, nu);
#endif
      }

    }
}
//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
//this one is used in GRID_SEARCH mode
//use individual frequencies in the parameter file
//so that spectral index can be calculated corresponding to the
//the two frequencies
void Evolution::results(Container *activeSh, std::vector<double> *nu,
                        double &tNow, std::vector<double> *fluxNu1,
                        std::vector<double> *fluxNu2)
{
  bool shell;
  double location;
  double ti;
  double lcIR=0.0; //lightcurve
  double lcRa=0.0;
  double shWidth;
  double innerRadius;
  double outerRadius;
  double bEneDens;
  double intEner;
  double specIndex;
  double logFluxRatio;
  double logNuRatio = log10(extnu2/extnu1);
  bool useResFile = extUseResFile;

   //std::cout<<setiosflags(std::ios::fixed)<<std::setprecision(8)<<tNow<<"in results"<<std::endl;
  std::string filename = extResultsFile;
  std::string lcFilename = extLCFile;
  std::ofstream timeFile(lcFilename.c_str(), std::ofstream::app);
  //timeFile.precision(8);

  if(timeFile.is_open())
    {
      timeFile<<tNow;

      for(Container::const_iterator sh = activeSh->begin();
          sh != activeSh->end(); ++sh)
        {
          shell = (*sh)->getShellId();
          std::vector<double> iNu = (*sh)->iNuVals();

          if(!shell)
            {
              lcRa += iNu[0]/1.e-26;
              lcIR += iNu[1]/1.e-26;
            }
        }

      if(lcIR == 0.0 || lcRa == 0.0)
        {
          logFluxRatio = 0.0;
        }
      else
        {
          logFluxRatio = log10(lcIR/lcRa);
        }

      specIndex = logFluxRatio/logNuRatio;
      fluxNu1->push_back(lcRa);
      fluxNu2->push_back(lcIR);

      timeFile<<"\t"<<lcRa<<"\t"<<lcIR<<"\t"<<specIndex<<"\n";

      lcIR = 0.0;
      lcRa = 0.0;

    }
  else
    {
      std::cerr<<"Error! Unable to open lightcurve file"<<std::endl;
    }
  timeFile.close();

if(useResFile)
    {
      std::ofstream resultsFile(filename.c_str(), std::ofstream::app);
      if (resultsFile.is_open())
        {

          for(Container::const_iterator sh = activeSh->begin();
              sh != activeSh->end(); ++sh)
          {
            shell = (*sh)->getShellId();
            location = (*sh)->getLocation();
            ti = (*sh)->getTimeOfInjection();
            innerRadius = (*sh)->getInnerRadius();
            outerRadius = (*sh)->getOuterRadius();
            shWidth = outerRadius - innerRadius;//(*sh)->getShellWidth();
            bEneDens = (*sh)->getBEneDens();
            intEner = (*sh)->getInternalEnergy();


            std::vector<double> iNu = (*sh)->iNuVals();

            if(location >= 0.0 && location <= 1.e20)
              {
                resultsFile<<"\n"<<location;
                if(!shell)
                  {
                    for(std::vector<double>::const_iterator it = iNu.begin(),
                          it3= nu->begin(); it != iNu.end(); ++it, ++it3)
                      {
                        resultsFile<<"\t"<<*it/1.e-26;
                      }
                    resultsFile<<"\t"<<shWidth;
                  }
                else
                  {
                    for(std::vector<double>::const_iterator it2 = nu->begin();
                        it2 != nu->end(); ++it2)
                      {
                        resultsFile<<"\t"<<0.0;

                      }
                    resultsFile<<"\t"<<shWidth;
                  }
              }
          }
          resultsFile<<"\n";
        }
      else
        {
          std::cout<<"Error! Unable to open results file"<<std::endl;
        }
      resultsFile.close();
    }
 }


//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
void Evolution::freqRange(std::vector<double> *nu)
{
  bool indivFreq = extindivFreq;
  double nu1 = extnu1;
  double nu2 = extnu2;

  if(!indivFreq)
    {
      double min = extnuMin;
      double max = extnuMax;

      int nuPoints = extnuPoints;

      double abscissa = min;
      double r = pow(max / min, 1.0 / nuPoints);
      double r2 = sqrt(r);
      r -= 1.0;

      //writing the frequency range file
      //useful when sampling many frequencies

      std::ofstream nuFile("nuRange.dat");
      std::cout<<"Log binning for frequency(Hz); written to nuRange.dat "
               <<std::endl;
      for (int i = 1; i <= nuPoints; ++i)
        {
          nu->push_back(r2 * abscissa); // Calc log-space mid-point for grid
          abscissa += (r * abscissa);  // Add grid width to abscissa to find
          //std::cout<<r2<<"\t"<<(abscissa)<<std::endl;
          nuFile<<r2 * abscissa<<std::endl;

        }
    }
  else
    {
      nu->push_back(nu1);
      nu->push_back(nu2);
    }


}


//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
void Evolution::avgLC(std::vector<double> *fNu1, std::vector<double> *fNu2)
{
  std::ofstream avgFile("avgSpecInd.dat", std::ofstream::app);
  std::ofstream paramFile("paramFile.dat", std::ofstream::app);

  std::vector<double>::const_iterator it, it2;
  int count = 1;
  double nu1Tot = 0.0;
  double nu2Tot = 0.0;
  double avgNu1 = 0.0;
  double avgNu2 = 0.0;
  double avgSPindex = 0.0;
  double logNuRatio = log10(extnu2/extnu1);

  double param1 = extJetLum;
  double param2 = extJetOpenAngle;
  double param3 = extEjTimeGap;
  double param4 = extGamMax;
  double param5 = extShWidth;
  double param6 = extShThermalEne;
  double param7 = extShIntKine;
  double param8 = extShIntMag;

  if(avgFile.is_open() && paramFile.is_open())
    {
      for (it = fNu1->begin(), it2 = fNu2->begin();
           it != fNu1->end(); ++it, ++it2)
        {
          count += 1;
          nu1Tot += *it;
          nu2Tot += *it2;
        }

      avgNu1 = nu1Tot/count;
      avgNu2 = nu2Tot/count;

      avgSPindex = (log10(nu2Tot/nu1Tot))/logNuRatio;

      avgFile<<avgNu1<<"\t"<<avgNu2<<"\t"
             <<avgSPindex<<std::endl;

      paramFile<<param1<<"\t"<<param2<<"\t"<<param3<<"\t"
               <<param4<<"\t"<<param5<<"\t"<<param6<<"\t"
               <<param7<<"\t"<<param8<<std::endl;

      count = 1;


    }
  else
    {
      std::cerr<<"Error! Unable to open file avgLC.dat"<<std::endl;
    }

  avgFile.close();

}
//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
void Evolution::finalTStepTau(Container *activeSh, std::vector<double> *nu,
                        double &tNow)
{
  bool shell;
  double location;
  double ti;
  double shWidth;
  double innerRadius;
  double outerRadius;
  double bEneDens;
  double intEner;

  std::ofstream finalTStauF("finalTStau.dat", std::ofstream::app);

  if (finalTStauF.is_open())
    {

      for(Container::const_iterator sh = activeSh->begin();
          sh != activeSh->end(); ++sh)
        {
          shell = (*sh)->getShellId();
          location = (*sh)->getLocation();
          ti = (*sh)->getTimeOfInjection();
          innerRadius = (*sh)->getInnerRadius();
          outerRadius = (*sh)->getOuterRadius();
          shWidth = outerRadius - innerRadius;//(*sh)->getShellWidth();
          bEneDens = (*sh)->getBEneDens();
          intEner = (*sh)->getInternalEnergy();


          std::vector<double> tauNu = (*sh)->tauVals();

          if(location >= 0.0 && location <= 1.e20)
            {
              finalTStauF<<std::setiosflags(std::ios::scientific)<<std::setprecision(8)<<location;
              if(!shell)
                {
                  for(std::vector<double>::const_iterator it3= nu->begin(), it4 = tauNu.begin();
                      it3 != nu->end();
                      ++it3, ++it4)
                    {
                      finalTStauF<<"\t"<<*it4;
                    }
                }
              else
                {
                  for(std::vector<double>::const_iterator it2 = nu->begin();
                      it2 != nu->end(); ++it2)
                    {
                      finalTStauF<<"\t"<<0.0;
                    }
                }
              finalTStauF<<"\n";
            }
        }
      finalTStauF<<"\n";

    }
  else
    {
      std::cout<<"Error! Unable to open finalTStau file"<<std::endl;
        }
  finalTStauF.close();


}
//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
void Evolution::finalTStep(Container *activeSh, std::vector<double> *nu,
                        double &tNow)
{
  bool shell;
  double location;
  double ti;
  double shWidth;
  double innerRadius;
  double outerRadius;
  double bEneDens;
  double intEner;

  std::ofstream finalTSFile("finalTS.dat", std::ofstream::app);

  if (finalTSFile.is_open())
    {

      for(Container::const_iterator sh = activeSh->begin();
          sh != activeSh->end(); ++sh)
        {
          shell = (*sh)->getShellId();
          location = (*sh)->getLocation();
          ti = (*sh)->getTimeOfInjection();
          innerRadius = (*sh)->getInnerRadius();
          outerRadius = (*sh)->getOuterRadius();
          shWidth = outerRadius - innerRadius;//(*sh)->getShellWidth();
          bEneDens = (*sh)->getBEneDens();
          intEner = (*sh)->getInternalEnergy();


          std::vector<double> iNu = (*sh)->iNuVals();

          if(location >= 0.0 && location <= 1.e20)
            {
              finalTSFile<<std::setiosflags(std::ios::scientific)<<std::setprecision(8)<<location;
              if(!shell)
                {
                  for(std::vector<double>::const_iterator it = iNu.begin(),
                        it3= nu->begin(); it != iNu.end(); ++it, ++it3)
                    {
                      finalTSFile<<"\t"<<*it/1.e-26;
                    }
                  finalTSFile<<"\t"<<shWidth;
                }
              else
                {
                  for(std::vector<double>::const_iterator it2 = nu->begin();
                      it2 != nu->end(); ++it2)
                    {
                      finalTSFile<<"\t"<<0.0;

                    }
                  finalTSFile<<"\t"<<shWidth;
                }
              finalTSFile<<"\n";
            }
        }
      finalTSFile<<"\n";
    }
  else
    {
      std::cout<<"Error! Unable to open results file"<<std::endl;
        }
  finalTSFile.close();

}
//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
bool Evolution::fileExists(const std::string& fileName)
{
  std::fstream fin;
  fin.open(fileName.c_str(),std::ios::in);
  if( fin.is_open() )
    {
      fin.close();
      return true;
    }
  else
    {
    fin.close();
    return false;
    }

}
//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
//this writes a file that contains every time step with
//with each shell location and optical depths. The optical depths
//will correspond to the frequencies given in nuRange.dat
void Evolution::tauWrite(Container *activeSh, std::vector<double> *nu)
{
  bool shell;
  double location;
  double ti;

  std::ofstream tauFile("tauResFile.dat", std::ofstream::app);

  if (tauFile.is_open())
    {

      for(Container::const_iterator sh = activeSh->begin();
          sh != activeSh->end(); ++sh)
        {
          shell = (*sh)->getShellId();
          location = (*sh)->getLocation();
          ti = (*sh)->getTimeOfInjection();

          std::vector<double> tauNu = (*sh)->tauVals();

          if(location >= 0.0 && location <= 1.e20)
            {
              tauFile<<std::setiosflags(std::ios::scientific)
                     <<std::setprecision(8)
                     <<location;
              if(!shell)
                {
                  for(std::vector<double>::const_iterator it3= nu->begin(),
                        it4 = tauNu.begin();
                      it3 != nu->end();
                      ++it3, ++it4)
                    {
                      tauFile<<"\t"<<*it4;
                    }
                }
              else
                {
                  for(std::vector<double>::const_iterator it2 = nu->begin();
                      it2 != nu->end(); ++it2)
                    {
                      tauFile<<"\t"<<0.0;
                    }
                }
              tauFile<<"\n";
            }
        }
      tauFile<<"\n";

    }
  else
    {
      std::cout<<"Error! Unable to open tauRes file"<<std::endl;
        }
  tauFile.close();
}
//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
void Evolution::splitLightcurve(Container *activeSh, double &tNow)
{
  std::string one="1", two="2", three="3", four="4", five="5", six="6";
  std::string lcFilename = extLCFile;
  std::string lcFilename1(lcFilename);
  std::string lcFilename2(lcFilename);
  std::string lcFilename3(lcFilename);
  std::string lcFilename4(lcFilename);
  std::string lcFilename5(lcFilename);
  std::string lcFilename6(lcFilename);
  int ss=lcFilename.size();
  int ss2=ss-4;
  lcFilename1.insert(ss2,one);
  lcFilename2.insert(ss2,two);
  lcFilename3.insert(ss2,three);
  lcFilename4.insert(ss2,four);
  lcFilename5.insert(ss2,five);
  lcFilename6.insert(ss2,six);

  std::ofstream timeFile(lcFilename.c_str(), std::ofstream::app);
  std::ofstream timeFile1(lcFilename1.c_str(), std::ofstream::app);
  std::ofstream timeFile2(lcFilename2.c_str(), std::ofstream::app);
  std::ofstream timeFile3(lcFilename3.c_str(), std::ofstream::app);
  std::ofstream timeFile4(lcFilename4.c_str(), std::ofstream::app);
  std::ofstream timeFile5(lcFilename5.c_str(), std::ofstream::app);
  std::ofstream timeFile6(lcFilename6.c_str(), std::ofstream::app);

  std::vector<double> nuVals(nu->size(), 0.0);
  std::vector<double> nuVals1(nu->size(), 0.0);
  std::vector<double> nuVals2(nu->size(), 0.0);
  std::vector<double> nuVals3(nu->size(), 0.0);
  std::vector<double> nuVals4(nu->size(), 0.0);
  std::vector<double> nuVals5(nu->size(), 0.0);
  std::vector<double> nuVals6(nu->size(), 0.0);

  std::vector<double>::iterator nuIt, nuIt1, nuIt2, nuIt3, nuIt4, nuIt5, nuIt6;
  std::vector<double>::const_iterator nuValIt, nuItConst, nuItConst1,
    nuItConst2, nuItConst3, nuItConst4, nuItConst5, nuItConst6;
  double up1=extUp1, up2=extUp2, up3=extUp3, up4=extUp4, up5=extUp5;

  double location;
  bool shell;

  if(timeFile.is_open())
    {
      timeFile<<std::setiosflags(std::ios::scientific)<<std::setprecision(8)<<tNow;
      timeFile1<<std::setiosflags(std::ios::scientific)<<std::setprecision(8)<<tNow;
      timeFile2<<std::setiosflags(std::ios::scientific)<<std::setprecision(8)<<tNow;
      timeFile3<<std::setiosflags(std::ios::scientific)<<std::setprecision(8)<<tNow;
      timeFile4<<std::setiosflags(std::ios::scientific)<<std::setprecision(8)<<tNow;
      timeFile5<<std::setiosflags(std::ios::scientific)<<std::setprecision(8)<<tNow;
      timeFile6<<std::setiosflags(std::ios::scientific)<<std::setprecision(8)<<tNow;

      for(Container::const_iterator sh = activeSh->begin();
          sh != activeSh->end(); ++sh)
        {
          location = (*sh)->getLocation();
          shell = (*sh)->getShellId();
          std::vector<double> iNu = (*sh)->iNuVals();

          if(!shell)
            {//total lightcurve
              for(nuValIt = iNu.begin(), nuIt = nuVals.begin();
                  nuValIt != iNu.end(); ++nuIt, ++nuValIt)
                {
                  *nuIt += *nuValIt/1.e-26;
                }
              //lightcurve1
              if(location <= up1)
                {
                  for(nuValIt = iNu.begin(), nuIt1 = nuVals1.begin();
                      nuValIt != iNu.end(); ++nuIt1, ++nuValIt)
                    {
                      *nuIt1 += *nuValIt/1.e-26;
                    }
                }//lightcurve2
              else if(location > up1 && location <= up2)
                {
                  for(nuValIt = iNu.begin(), nuIt2 = nuVals2.begin();
                      nuValIt != iNu.end(); ++nuIt2, ++nuValIt)
                    {
                      *nuIt2 += *nuValIt/1.e-26;
                    }
                }//lightcurve3
              else if(location > up2 && location <= up3)
                {
                  for(nuValIt = iNu.begin(), nuIt3 = nuVals3.begin();
                      nuValIt != iNu.end(); ++nuIt3, ++nuValIt)
                    {
                      *nuIt3 += *nuValIt/1.e-26;
                    }
                }//lightcurve4
              else if(location > up3 && location <= up4)
                {
                  for(nuValIt = iNu.begin(), nuIt4 = nuVals4.begin();
                      nuValIt != iNu.end(); ++nuIt4, ++nuValIt)
                    {
                      *nuIt4 += *nuValIt/1.e-26;
                    }
                }//lightcurve5
              else if(location > up4 && location <= up5)
                {
                  for(nuValIt = iNu.begin(), nuIt5 = nuVals5.begin();
                      nuValIt != iNu.end(); ++nuIt5, ++nuValIt)
                    {
                      *nuIt5 += *nuValIt/1.e-26;
                    }
                }//lightcurve6
               else if(location > up5)
                 {
                  for(nuValIt = iNu.begin(), nuIt6 = nuVals6.begin();
                      nuValIt != iNu.end(); ++nuIt6, ++nuValIt)
                    {
                      *nuIt6 += *nuValIt/1.e-26;
                    }
                 }
            }
        }

      for(nuItConst = nuVals.begin(), nuItConst1 = nuVals1.begin(),
            nuItConst2 = nuVals2.begin(), nuItConst3 = nuVals3.begin(),
            nuItConst4 = nuVals4.begin(),
            nuItConst5 = nuVals5.begin(), nuItConst6 = nuVals6.begin();
          nuItConst != nuVals.end(); ++nuItConst, ++nuItConst1, ++nuItConst2,
            ++nuItConst3, ++nuItConst4, ++nuItConst5, ++nuItConst6)
        {
          timeFile<<"\t"<<*nuItConst;
          timeFile1<<"\t"<<*nuItConst1;
          timeFile2<<"\t"<<*nuItConst2;
          timeFile3<<"\t"<<*nuItConst3;
          timeFile4<<"\t"<<*nuItConst4;
          timeFile5<<"\t"<<*nuItConst5;
          timeFile6<<"\t"<<*nuItConst6;
        }

      timeFile<<"\n";
      timeFile1<<"\n";
      timeFile2<<"\n";
      timeFile3<<"\n";
      timeFile4<<"\n";
      timeFile5<<"\n";
      timeFile6<<"\n";
      nuVals.clear();
      nuVals1.clear();
      nuVals2.clear();
      nuVals3.clear();
      nuVals4.clear();
      nuVals5.clear();
      nuVals6.clear();
    }
  else
    {
      std::cerr<<"Error! Unable to open lightcurve file"<<std::endl;
    }
  timeFile.close();
  timeFile1.close();
  timeFile2.close();
  timeFile3.close();
  timeFile4.close();
  timeFile5.close();
  timeFile6.close();

}

//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
void Evolution::compTable(std::vector<double> *nu, std::vector<double> *g)
{
  std::vector<double>::iterator nuIter_1, nuIter_0, gIter;

  try
    {
      // Create the xF_c(x) array as size_nu*size_nu*size_gamma
      comptonVect = new std::vector<double>((*nu).size()*(*nu).size()*(*g).size(),
				       0.0);
      for(gIter = g->begin(); gIter != g->end(); ++gIter)
	{
	  for(nuIter_1 = nu->begin(); nuIter_1 != nu->end(); ++nuIter_1)
	    {
	      for(nuIter_0 = nu->begin(); nuIter_0 != nu->end(); ++nuIter_0)
		{



		}
	    }
	}

    }
  catch(std::bad_alloc)
    {
      std::cerr<<"Memory error in compTable while creating the table";
    }
}
//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
void Evolution::gammaRange(std::vector<double> *egamma)
{
  //egamma = new std::vector<double>;
  double min = exteGammaMin;
  double max = exteGammaMax;

  int gammaPoints = extgammaPoints;

  double abscissa = min;
  double r = pow(max / min, 1.0 / gammaPoints);
  double r2 = sqrt(r);
  r -= 1.0;

  //writing the frequency range file
  //useful when sampling many frequencies

  std::ofstream gammaFile("gammaRange.dat");
  std::cout<<"Log binning for electrons(Lorentz factor)"
	   <<std::endl;
  for (int i = 1; i <= gammaPoints; ++i)
    {
      egamma->push_back(r2 * abscissa); // Calc log-space mid-point for grid
      abscissa += (r * abscissa);  // Add grid width to abscissa to find
      //std::cout<<r2<<"\t"<<(abscissa)<<std::endl;
      gammaFile<<r2 * abscissa<<std::endl;

    }

}
