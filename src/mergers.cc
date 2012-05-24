//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo....//
//                                                                    //
// Internal shock model for relativistic astrophysical jets           //
//                                                                    //
// Author: Omar Jamil                                                 //
// Contact: o.jamil@phys.soton.ac.uk                                  //
//          University of Southampton                                 //
//                                                                    //
//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo....//

//#include <functional>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <fstream>
#include "mergers.hh"
#include "funcobj.hh"
#include "physcon.hh"
#include "paraccess.hh"

//a small macro for square function
#define sqr(a) ((a)*(a))

//Merger treatment shown in Spada '01 and Panaitescu '99 followed
//with some modifications.

//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....

Mergers::Mergers(Shells *inner_, Shells *outer_) :
  innerShell(inner_), outerShell(outer_)
{}

//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....

Mergers::~Mergers()
{
 //  delete innerShell;
//   delete outerShell;
  
}
//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
void Mergers::doMerger(Container *active, double &tNow)
{
  try
    {
      
      
      Container::iterator inner, outer, locator, newShell;
      Container::iterator innerMinus, outerPlus;
      
      double adiabaticInd = 4./3.;
      double t_inner = 0.0, t_outer=0.0;
     
      t_inner = this->getInnerShell();
      t_outer = this->getOuterShell();
      
      locator = std::find_if(active->begin(), active->end(), 
                             FindShellWithDouble(t_inner));
      
      inner = std::find_if(active->begin(), active->end(), 
                           FindShellWithDouble(t_inner));
      outer = std::find_if(active->begin(), active->end(),
                           FindShellWithDouble(t_outer));
      
  
      
      
     //  (*inner)->setDimensions(tNow);
//       (*outer)->setDimensions(tNow);
      bool insh = (*inner)->getShellId();
      bool oush = (*outer)->getShellId();
      
      double innIntEne = (*inner)->getInternalEnergy();
      double outIntEne = (*outer)->getInternalEnergy();
      double innGamma = (*inner)->getShellGamma();
      double outGamma = (*outer)->getShellGamma();
      double innBeta = (*inner)->getShellBeta();
      double outBeta = (*outer)->getShellBeta();
      double innMass = (*inner)->getShellMass();
      double outMass = (*outer)->getShellMass();
      int innMergCount = (*inner)->getMergeCount();
      int outMergCount = (*outer)->getMergeCount();
      int totMergCount = innMergCount + outMergCount;
      double innOuterRadius = (*inner)->getOuterRadius();
      double innInnerRadius = (*inner)->getInnerRadius();
      double outInnerRadius = (*outer)->getInnerRadius();
      double outOuterRadius = (*outer)->getOuterRadius();
      double innMu = (innMass + innIntEne/sqr(physcon.c));
      double outMu = (outMass + outIntEne/sqr(physcon.c));
      
      double innWidth = (*inner)->getShellWidth();
      double outWidth = (*outer)->getShellWidth();
      
      
      //the following condition can cause a seg faul
      //allow the code to exit normally.
      if(innGamma == outGamma)
        {
          if(innIntEne == 0.0 && outIntEne == 0.0)
            {
              std::cerr<<"** Problem in doMergers::Mergers.cc"<<"\n"
                       <<"** Cannot have two shells with equal BLFs and zero int Ene colliding"<<"\n"
                       <<"** Segmentation fault imminent so I shall exit normally" <<"\n"
                       <<"Solution: E_frac > 0.0, or EThermal_frac = 0.0, or have shells with different BLFs"<<"\n"<<"\n"
                       <<"exit(1) called"<<"\n";
              exit(1);
            }
          else
            {
              std::cerr<<"**Same BLF shells merging -- not ideal and stable"<<"\n"
                       <<"There will still be a compression effect on the new shell"<<"\n"
                       <<"further investigation necessary to calculate energetics"<<"\n";
              
            }
          
        }
      
      //on rare occasions they are out by a fraction!
      //down to variations because of the large numbers and accuracy
      //possible. Within reasonable error bounds
      if(innOuterRadius > outInnerRadius)
        {
          innOuterRadius = outInnerRadius;
        }
      std::cout<<std::setiosflags(std::ios::scientific)<<std::setprecision(5);
      
      //The merged shell lorentz factor
     //  double mergedGamma = sqrt(((innMu * innGamma) + (outMu * outGamma))/
//                                 ((innMu / innGamma) + (outMu / outGamma)));
      
      // double mergedBeta = sqrt(1. - (1. / sqr(mergedGamma)));
      
      double c = physcon.c;
      double c2 = sqr(physcon.c);
      double beta_m2 = sqr(((innGamma*innMass*innBeta) + (outGamma*outMass*outBeta)) / 
                        ((innGamma*innMass + outGamma*outMass)));
      double mergedBeta = sqrt(beta_m2);
      double gamma_m = (1./sqrt(1. - beta_m2));

      double mergedGamma = gamma_m;
      
      // std::cout<<beta_m2<<std::endl;
      
//       std::cout<<innGamma<<"\t"<<outGamma<<"\t"<<mergedGamma<<"\t"
//                <<gamma_m<<std::endl;
      
      //shocked material gamma in the unshocked inner shell frame
      //from Panaitescu et al. 1999
      double gamShoInnT = mergedGamma * innGamma * (1.-mergedBeta*innBeta);
      double betaShoInnT = sqrt(1. - (1. / sqr(gamShoInnT)));
      
      double betaRS =  ((gamShoInnT - 1.)*(adiabaticInd * gamShoInnT + 1.))/
        (betaShoInnT * gamShoInnT*(adiabaticInd*(gamShoInnT - 1.) + 1.));
  
      //slow energization; shock crossing time for reverse shock
      double tCrossR = innWidth / (betaRS * physcon.c);
      
      double gamShoInn = 1./sqrt(1.-(betaRS * betaRS));
      
      
      double gamShoOutT = mergedGamma * outGamma * (1.-(mergedBeta*outBeta));
      
      double betaShoOutT = sqrt(1. - (1. / sqr(gamShoOutT)));
     

      double betaFS =  ((gamShoOutT - 1.)*(adiabaticInd * gamShoOutT + 1.))/
        (betaShoOutT * gamShoOutT*(adiabaticInd*(gamShoOutT - 1.) + 1.));
      
      //slow energization; shock crossing time for forward shock
      double tCrossF = outWidth / (betaFS * physcon.c);
      
      //slow energization; total crossing time for the shocks
      double tCrossT = tCrossR + tCrossF;
      
      double gamShoOut =  1./sqrt(1.-(betaFS * betaFS)); 
            
      double inMergedShellWidth = inNewShellWidth(adiabaticInd, 
                                                  gamShoInn, innGamma,
                                                  mergedGamma, 
                                                  innWidth);
      
      double ouMergedShellWidth = ouNewShellWidth(adiabaticInd, 
                                                  gamShoOut, outGamma,
                                                  mergedGamma, 
                                                  outWidth);
      
      double mergedTotalWidth = inMergedShellWidth + ouMergedShellWidth;
      
      //Merged Shell inner Radius
      double mergeInnerRadius = innOuterRadius-inMergedShellWidth;
      double mergeOuterRadius = innOuterRadius+ouMergedShellWidth;
            
      //std::cout<<inMergedShellWidth<<"\t"<<innWidth<<"\t"<<"inwidths"<<std::endl;
      //std::cout<<ouMergedShellWidth<<"\t"<<outWidth<<"\t"<<"ouwidths"<<std::endl;
      //std::cout<<innOuterRadius<<"\t"<<outInnerRadius<<"\t"<<"inn out"<<std::endl;
      //std::cout<<innInnerRadius<<"\t"<<outOuterRadius<<"\t"<<"innIn outOu"<<std::endl;
      //std::cout<<mergeInnerRadius<<"\t"<<mergeOuterRadius<<"\t"<<"Radii"<<std::endl;
      //std::cout<<"Diff"<<"\t"<<totTemp - mergedTotalWidth<<std::endl;
      
      
      double shCentre = 0.5 * (mergeInnerRadius + mergeOuterRadius);
      //Based on Merger Gamma the merged shell Time of Inj
      double mergedTimeInj = tNow - (shCentre/
                                     (mergedBeta*physcon.c));
      
     //  std::cout<<t_inner<<"\t"<<t_outer<<"\t"<<mergedTimeInj
//                <<"\t"<<shCentre<<std::endl;
      
      //slow energization; time of full of energization
      double tEner = tNow + tCrossT;
      
      double mergedMass = innMass + outMass;
      double mergedIntEne = mergedShellIntEner(innIntEne, outIntEne,
                                               innMu, outMu,
                                               innGamma, outGamma,
                                               mergedGamma);
      bool slowEnergization = extSlowEnerg;
      double intEner = 0.0;
  
      if(slowEnergization)
        {
          //for slowed energization; merged shell inserted with zero
          //internal energy and then slowly increased
          active->insert(locator, new Shells(mergedTimeInj, mergedMass,
                                         mergedGamma, intEner,
                                             mergedTotalWidth));
        }
      else
        {
          active->insert(locator, new Shells(mergedTimeInj, mergedMass,
                                             mergedGamma, mergedIntEne,
                                             mergedTotalWidth));
        }
      
      
      //Set the upper and lower radius of the newly inserted shell
      newShell =std::find_if(active->begin(), active->end(), 
                             FindShellWithDouble(mergedTimeInj));
            
      (*newShell)->setLocation(shCentre);
 
//****TEMPORARY??****   
  bool sp = extSplitLC;
  if(sp)
    {
      mergeLocCount(shCentre);
    }
//****??TEMPORARY****
    
      (*newShell)->setLowerAndUpperX(mergeInnerRadius, 
                                     mergeOuterRadius);
      
      //setting the shell id false i.e. a merged shell
      (*newShell)->setIdFalse();
      //increasing the merger count
      (*newShell)->setMergeCount(totMergCount);
      //setting the merged shell volume
      (*newShell)->setInitialMergedVA();
      
      //calculating the powerlaw normalization for the merged shell
      
      (*newShell)->powerlawNorm();

      //(*newShell)->setExpansionBeta();
      (*newShell)->initialBEnergyDensity();
      (*newShell)->setPrevTime(tNow);
      
      if(slowEnergization)
        {
          //slow energization related
          (*newShell)->setTEnergize(tEner);
          (*newShell)->setIntEavail(mergedIntEne);
        }
      
      //removing the merged shells from the list
      delete *inner;
      delete *outer;
      active->erase(inner);
      active->erase(outer); 
      
      
      //std::cout<<(*locator)->getTimeOfInjection()<<"\n";
    }
  catch(std::bad_alloc)
    {
      std::cerr<<"Memory exhauseted in doMerger;Merger.cc\n";
    }  
}
//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
double Mergers::inNewShellWidth(double &adia, double &gShIn, 
                                double &gaIn, double &gaMer, 
                                double &widIn)
{
  //from Spada '01 paper but with extra Lorentz factors to
  //calculate in the lab frame
  //the number density spada eq 10
  double rhoIn =  (gaMer/gaIn)*
    (((adia*gShIn) + 1.)/(adia - 1.));
  
  //double rhoIn = (((adia*gShIn) + 1.)/(adia - 1.));
  //Spada equation 13
  double deltaIn = widIn/rhoIn;
  
  double InDeltaNewShell = deltaIn;
  //double InDeltaNewShell = widIn;
  
  return InDeltaNewShell;
  
  
}

//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
double Mergers::ouNewShellWidth(double &adia, double &gShOu,  
                                double &gaOu, double &gaMer,
                                double &widOu)
{
  //from Spada '01 paper but with extra Lorentz factors to
  //calculate in the lab frame
  //the number density spada eq 10
  double rhoOu = (gaMer/gaOu)*
    (((adia*gShOu) + 1.)/(adia - 1.));
   
  //double rhoOu = (((adia*gShOu) + 1.)/(adia - 1.));
  //Spada equation 13
  double deltaOu = widOu/rhoOu;
  
  double OuDeltaNewShell = deltaOu;
  //double OuDeltaNewShell = widOu;
  
  return OuDeltaNewShell;
  
  
}
//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
double Mergers::mergedShellIntEner(double &InEta, double &OuEta,
                                   double &InMu, double &OuMu,
                                   double &InGam, double &OuGam,
                                   double &MerGam)
{

  //this is calculated in the lab frame.

 //  double mergeInternalEnergy = InEta + OuEta + 
//     (InMu*sqr(physcon.c)*(InGam - MerGam)) +
//     (OuMu*sqr(physcon.c)*(OuGam - MerGam));

   double mergeInternalEnergy = MerGam*InEta + MerGam*OuEta + 
    (InMu*sqr(physcon.c)*(InGam - MerGam)) +
    (OuMu*sqr(physcon.c)*(OuGam - MerGam));
  
  // std::cout<<InEta<<"\t"<<OuEta<<"\t"<<InMu<<"\t"<<OuMu<<"\t"
//            <<InGam-MerGam<<"\t"<<OuGam-MerGam<<"\t"<<MerGam<<"\t"
//            <<mergeInternalEnergy<<"\t"<<"meR"<<std::endl;
  
  if(mergeInternalEnergy < 0.0)
    {
      mergeInternalEnergy = 0.0;
    }
    
  return mergeInternalEnergy;
  
}
//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
void Mergers::mergeLocCount(double &location)
{
  std::ofstream locCount("locCount.dat", std::ofstream::app);
  
  double up1=extUp1, up2=extUp2, up3=extUp3, up4=extUp4, up5=extUp5;
  if (locCount.is_open()) 
    { 
  if(location <= up1)
    {
      locCount <<"1"<<"\t"<<"0"<<"\t"<<"0"<<"\t"<<"0"
               <<"\t"<<"0"<<"\t"<<"0"<<std::endl;
    }
  else if(location > up1 && location <= up2)
    {
      locCount <<"0"<<"\t"<<"1"<<"\t"<<"0"<<"\t"<<"0"
               <<"\t"<<"0"<<"\t"<<"0"<<std::endl;
    }//lightcurve3
  else if(location > up2 && location <= up3)
    {
     locCount <<"0"<<"\t"<<"0"<<"\t"<<"1"<<"\t"<<"0"
               <<"\t"<<"0"<<"\t"<<"0"<<std::endl;
    } 
    //lightcurve4
  else if(location > up3 && location <= up4)
    {
      locCount <<"0"<<"\t"<<"0"<<"\t"<<"0"<<"\t"<<"1"
               <<"\t"<<"0"<<"\t"<<"0"<<std::endl;
    }
    //lightcurve5
  else if(location > up4 && location <= up5)
    {
      locCount <<"0"<<"\t"<<"0"<<"\t"<<"0"<<"\t"<<"0"
               <<"\t"<<"1"<<"\t"<<"0"<<std::endl;
    }
    //lightcurve6
  else if(location > up5)
    {
      locCount <<"0"<<"\t"<<"0"<<"\t"<<"0"<<"\t"<<"0"
               <<"\t"<<"0"<<"\t"<<"1"<<std::endl;
    }
 }
  else 
    {
      std::cout<<"Error! Location count file error"<<std::endl;
        }
  locCount.close();   
}
//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
