//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo....//
//                                                                    //
// Internal shock model for relativistic astrophysical jets           //
//                                                                    //
// Author: Omar Jamil                                                 //
// Contact: o.jamil@phys.soton.ac.uk                                  //
//          University of Southampton                                 //
//                                                                    //
//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo....//

#include "parameters.hh"
#include "numerical.hh"
#include "physcon.hh"

//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
Parameters::Parameters(std::string fn)
{
  try 
    {
      loadParametersFile(fn);
    } 
  catch (const std::exception &e) 
    {
      throw;
    }
    
}

//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
Parameters::~Parameters()
{}

//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
void Parameters::loadParametersFile(std::string filename)
  throw (const std::exception &)
{
  std::ifstream parameterFile(filename.c_str());
  if (!parameterFile.good())
    throw std::runtime_error("Parameter file failed to load");

  std::string label, dummy, dummy2;
  std::cout<<"Reading the parameters file...."<<"\n";
  
  
  while(!parameterFile.eof())
    {
      parameterFile >> label;
      
      if(label == "load_prev_shells")
        {
          parameterFile >> dummy >> dummy2;
          if(dummy2 == "y")
            {
              loadPrevSim_ = 1;
            }
          else
            {
              loadPrevSim_ = 0;
            }
        }
      else if(label == "dump_active_at_end")
        {
          parameterFile >> dummy >> dummy2;
          if(dummy2 == "y")
            {
              dumpActive_ = 1;
            }
          else
            {
              dumpActive_ = 0;
            }
        }
      else if(label == "jet_luminosity")
        {
          parameterFile >> dummy >> jetLuminosity_;
        }
      else if(label == "jet_opening_angle")
        {
          parameterFile >> dummy >> openingAngle_;
        }
      else if(label == "jet_viewing_angle")
        {
          parameterFile >> dummy >> viewingAngle_;
        }
      else if(label == "source_distance")
        {
          parameterFile >> dummy >> sourceDistance_;
        }
      else if(label == "avg_ejection_gap")
        {
          parameterFile >> dummy >> ejectionGap_;
        }
      else if(label == "shell_inj_duration")
        {
          parameterFile >> dummy >> duration_;
        }
      else if(label == "BLF_min")
        {
          parameterFile >> dummy >> GammaMin_;
        }
      else if(label == "BLF_max")
        {
          parameterFile >> dummy >> GammaMax_;
        }
      else if(label == "shell_width_factor")
        {
          parameterFile >> dummy >> width_;
        }
      else if(label == "EThermal_frac")
        {
          parameterFile >> dummy >> thermalEne_;
        }
      else if(label == "EelecKin_frac")
        {
          parameterFile >> dummy >> EintKine_;
        }
      else if(label == "EMagnet_frac")
        {
          parameterFile >> dummy >> EMagnetic_;
        }
      else if(label == "powerlaw_index")
        {
          parameterFile >> dummy >> pLawIndex_;
        }
      else if(label == "e_gamma_min")
        {
          parameterFile >> dummy >> gMin_;
        }
      else if(label == "e_gamma_max")
        {
          parameterFile >> dummy >> gMax_;
        }
      else if(label == "nu_min")
        {
          parameterFile >> dummy >> nuMin_;
        }
      else if(label == "nu_max")
        {
          parameterFile >> dummy >> nuMax_;
        }
      else if(label == "nu_points")
        {
          parameterFile >> dummy >> nuPoints_;
        }
      else if(label == "gamma_points")
        {
          parameterFile >> dummy >> gammaPoints_;
        }
      else if(label == "individual_frequencies")
        {
          parameterFile >> dummy >> dummy2;
          if(dummy2 == "y")
            {
              individual_frequencies_ = 1;
            }
          else
            {
              individual_frequencies_ = 0;
            }
        }
      else if(label == "nu_1")
        {
          parameterFile >> dummy >> nu_1_;
        }
      else if(label == "nu_2")
        {
          parameterFile >> dummy >> nu_2_;
        }
      else if(label == "use_shell_file")
        {
          parameterFile >> dummy >> dummy2;
          if(dummy2 == "y")
            {
              useShellFile_ = 1;
            }
          else
            {
              useShellFile_ = 0;
            }
        }
      else if(label == "shell_file")
        {
          parameterFile >> dummy >> shellFile_;
        }
      else if(label == "use_data_file")
        {
          parameterFile >> dummy >> dummy2;
          if(dummy2 == "y")
            {
              useDataFile_ = 1;
            }
          else
            {
              useDataFile_ = 0;
            }
        }
      else if(label == "data_file")
        {
          parameterFile >> dummy >> dataFile_;
        }
      else if(label == "write_results_file")
        {
          parameterFile >> dummy >> dummy2;
          if(dummy2 == "y")
            {
              useResultsFile_ = 1;
            }
          else
            {
              useResultsFile_ = 0;
            }
        }
      else if(label == "results_file")
        {
          parameterFile >> dummy >> resultsFile_;
        }
      else if(label == "write_every")
        {
          parameterFile >> dummy >> everyTStep_;
        }
       else if(label == "final_time_step")
        {
          parameterFile >> dummy >> dummy2;
          if(dummy2 == "y")
            {
              finalTStep_ = 1;
            }
          else
            {
              finalTStep_ = 0;
            }
        }
      else if(label == "lightcurve_file")
        {
          parameterFile >> dummy >> lcFile_;
        }
      else if(label == "increase_time_resolution")
        {
          parameterFile >> dummy >> dummy2;
          if(dummy2 == "y")
            {
              increaseTimeStep_ = 1;
            }
          else
            {
              increaseTimeStep_ = 0;
            }
        }
      else if(label == "step_resolution")
        {
          parameterFile >> dummy >> stepResolution_;
        }
      else if(label == "total_run_duration")
        {
          parameterFile >> dummy >> totalDuration_;
        }
       else if(label == "in_vacuum_expansion")
        {
          parameterFile >> dummy >> dummy2;
          if(dummy2 == "y")
            {
              inVacuum_ = 1;
            }
          else
            {
              inVacuum_ = 0;
            }
        }
      else if(label == "inj_int_energy")
        {
          parameterFile >> dummy >> dummy2;
          if(dummy2 == "y")
            {
              injIntEne_ = 1;
            }
          else
            {
              injIntEne_ = 0;
            }
        }
      else if(label == "rel_mass_frac")
        {
          parameterFile >> dummy >> relMassFrac_;
        }
      else if(label == "slow_energization")
        {
          parameterFile >> dummy >> dummy2;
          if(dummy2 == "y")
            {
              slowEnergization_ = 1;
            }
          else
            {
              slowEnergization_ = 0;
            }
        }
      else if(label == "shock_location")
        {
          parameterFile >> dummy >> shocLoc_;
        }
    }

  //****TEMPORARY??****
  std::ifstream splitPar("splitLC.par");
  std::string lab, dum;
  while(!splitPar.eof())
    {
      splitPar >> lab;
      if(lab == "split?")
        {
          splitPar >> dum >> split;
        }
      if(lab == "upper_1")
        {
           splitPar >> dum >> up1;
        }
      if(lab == "upper_2")
        {
          splitPar >> dum >> up2;
        }
       if(lab == "upper_3")
        {
          splitPar >> dum >> up3;
        }
       if(lab == "upper_4")
        {
          splitPar >> dum >> up4;
        }
       if(lab == "upper_5")
        {
          splitPar >> dum >> up5;
        } 
             
    }
  //****??TEMPORARY****
}

//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
double Parameters::synchConstap()
{
  Numerical gFun;
  double p = pLawIndex_;
  double temp1 = p/4. + 19./12.;
  double temp2 = p/4. - 1./12.;
  double temp3 = p/4. + 5./4.;
  double temp4 = p/4. + 7./4.;
  
  double a_p = (sqrt(3.14159265359)/2.)*((gFun.gammaFunction(temp1)*
                                          gFun.gammaFunction(temp2)
                                          *gFun.gammaFunction(temp3))/
                                         ((p+1.)*gFun.gammaFunction(temp4)));
  return a_p;
  
}
//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
double Parameters::synchConstbp()
{
  Numerical gFun;
  double p = pLawIndex_;
  double temp1 = (3.*p+22.)/12.;
  double temp2 = (3.*p+2.)/12.;
  double temp3 = (p+6.)/4.;
  double temp4 = (p+8.)/4.;
  
  double b_p = (sqrt(3.14159265359)/8.)*((gFun.gammaFunction(temp1)*
                                          gFun.gammaFunction(temp2)
                                          *gFun.gammaFunction(temp3))/
                                         (gFun.gammaFunction(temp4)));
  return b_p;
  
}
//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
