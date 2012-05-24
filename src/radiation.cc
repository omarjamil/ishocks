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

#include "radiation.hh"
#include "paraccess.hh"
#include "physcon.hh"

//a small macro for square function
#define sqr(a) ((a)*(a))

//Ref: Longair, '02 ~pg 250 - 262
//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
Radiation::Radiation()
{}

//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
Radiation::~Radiation()
{}

//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....

void Radiation::synchSpectrum(std::vector<double> *nu, Container *actSh)
{

  std::vector<double>::const_iterator nuIter, nuIter2, gIter, eDistIter;
  double Jnu, Inu, shGamma, shTau;
  double absor, emiss;
  double intEn, kap;
  double shVol, shArea, shWidth, shOuterRad;
  bool shell;
  double shTi;
  double shDoppler;
  double nuEmit;
  double shBDens;
  double shGMax;
  double nuCrit;
  double nuGyro;
  double shMass;
  double distance = extSourceDistance;
  double expAbs;
  double shellLum0=0.0, totalInt=0.0, jetLum=0.0, Inu2=0.0;
  std::vector<double> *gamma, *eDist;
  
      
  for (Container::const_iterator sIt = actSh->begin();
       sIt != actSh->end(); ++sIt)
    {
      shell = (*sIt)->getShellId();
      shGamma = (*sIt)->getShellGamma();
      intEn = (*sIt)->getInternalEnergy();
      kap = (*sIt)->getPLawNorm();
      shVol = (*sIt)->getShellVolume();
      shArea = (*sIt)->getShellArea();
      shWidth = (*sIt)->getShellWidth();
      shOuterRad = (*sIt)->getShellOuRad();
      shTi = (*sIt)->getTimeOfInjection();
      shDoppler = (*sIt)->dopplerFactor();
      shBDens = (*sIt)->getBEneDens();
      shGMax = (*sIt)->getGammaMax();
      shMass = (*sIt)->getShellMass();
      nuGyro = (physcon.e * shBDens)/(2.*physcon.pi*physcon.m_e);
      nuCrit = sqr(shGMax) * nuGyro;
     
      
      if(!shell)
        {
          
          (*sIt)->nuClear();
          (*sIt)->tauClear(); 

          for(nuIter = nu->begin(); nuIter != nu->end(); ++nuIter)
            {
              //the frequency in shell rest frame
              nuEmit = (*nuIter)/shDoppler;
              //nuEmit = *nuIter;
              //nuEmit = (*nuIter);
              if(nuEmit >= nuCrit)
                {
                  Inu = 0.0;
                }
              else 
                {
                    
                  emiss = synchrotronEmiss(intEn, nuEmit, kap, shBDens);
                  absor = synchrotronAbsor(intEn, nuEmit, kap, shBDens);
              
                  //anything less is a numerical error
                  if(absor <= 1.e-99 || emiss <= 1.e-99) 
                    {
                      Inu = 0.0;
                    }
                  else
                    {
                  
                      expAbs = exp(-absor*shWidth);
              
                      if(expAbs <= 1.e-99)
                        {
                          expAbs = 0.0;
              
                        }
                      
                      //Doppler factor to observe in lab frame
                      Inu =  (1./(physcon.pi*sqr(distance)))*
                        pow(shDoppler,3)*shArea * 
                        ((emiss)/(absor * 4. * physcon.pi))
                        *(1. - exp(-(absor*shOuterRad)));
                      shTau = absor*shOuterRad;
                      //std::cout<<shTi<<"\t"<<*nuIter<<"\t"<<nuEmit<<"\t"<<intEn
                      //         <<"\t"<<shOuterRad<<"\t"<<shArea<<"\t"<<Inu/1.e-26<<std::endl;
                      Inu2 = shArea*((emiss)/(absor * 4. * physcon.pi))
                        *(1. - exp(-(absor*shOuterRad)));
                      shellLum0+=(*nuIter*Inu2);
                      //std::cout<<*nuIter<<"\t"<<Inu<<"\t"<<*nuIter*Inu<<std::endl;
                      
                      
                      emiss = 0.0;
                      absor = 0.0;
                      
                      if (Inu < 0.0)
                        {
                          Inu = 0.0;
                        }
                      
                            
                    }
                                            
                }
              
               (*sIt)->nuWrite(Inu);
               (*sIt)->tauWrite(shTau);

               Inu2 = 0.0;
               shTau = 0.0;
            }
          if(shellLum0 > (shGamma*shMass*9.e16))
            {
              std::cerr<<"****Shell "<<shTi<<" radiative Luminosity too high!"<<"\n"
                       <<"****Luminosity: (W)"<<"\t"<<shellLum0<<"\n"
                       <<"****Shell Relativistic Energ (J)"<<"\t"<<(shGamma*shMass*9.e16)<<"\n"
                       <<"****Hit radiative losses regime; reduce BLF range..."<<"\n"
                       <<"****rel_mass_frac, increase shell size"<<"\n"
                       <<"*****************************************************"<<std::endl;
            }
          
          jetLum+=shellLum0;
          totalInt+=intEn;
          shellLum0=0.0;
          
          
              
        } else 
        {
          continue;
        }
     
    }
  if(jetLum > extJetLum)
    {
      std::cerr<<"****Jet Luminosity too high!"<<"\n"
               <<"****Jet Luminosity: (W)"<<"\t"<<jetLum<<"\n"
               <<"****Total Int. Ene.(J)"<<"\t"<<totalInt<<"\n"
               <<"****Hit radiative losses regime; reduce BLF range..."<<"\n"
               <<"****rel_mass_frac, increase shell size"<<"\n"
               <<"*****************************************************"<<std::endl;
    }
  jetLum=0.0;
  totalInt=0.0;
  
}

//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
double Radiation::synchrotronEmiss(double &intEne, double &nu, 
                                   double &kappa, double &BDen)
{
  
  //N(E) electrons per m^-3 per mc^2 energy in GeV
  //from the calculation of powerlaw normalization
  
  
  double a_p = extSynchConstAp;
  double BField = BDen;
  double p = extPowLawInd;
  double temp1 = (p+1.)/2.;
  double temp2 = (p-1.)/2.;
  

//   double JNu = (2.344e-25)*a_p*pow(BField, temp1)
//     *kappa*pow((3.217e17/nu), temp2);

  //N(E) in SI units 
  double JNu = (2.344e-25)*a_p*pow(BField, temp1)
    *kappa*pow(1.253e37/nu, temp2);

//   double JNu = (pow(3,0.5)*BField*kappa)/(4.*physcon.pi*physcon.epsilon_0*
//                                           physcon.c*physcon.m_e)
//     *pow(((3*physcon.e*BField)/(2*physcon.pi*nu*pow(physcon.m_e,3)*
//                                 pow(physcon.c,4))),temp2)*a_p;

  
  
  //units W m^-3 Hz^-1
  return JNu; 

}
//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
double Radiation::synchrotronAbsor(double &intEne, double &nu, 
                                   double &kappa, double &BDen)
{
  //N(E) electrons per m^-3 per mc^2 energy in GeV
  //from the calculation of powerlaw normalization

  double b_p = extSynchConstBp;
  double BField = BDen;
  double p = extPowLawInd;
  double temp1 = (p+2.)/2.;
  double temp2 = (p+4.)/2.;
  
  //Longair pg 262 
 //  double ChiNu = 20.9*kappa*pow(BField,temp1)
//     *pow((5.67e9),p)*b_p*pow(nu,-temp2);

  //N(E) in SI units 
  double ChiNu = 3.354e-9*kappa*pow(BField,temp1)
    *pow((3.354e18),p)*b_p*pow(nu,-temp2);

//   double ChiNu = (pow(3,0.5)*pow(physcon.e,3)*physcon.c)/(8.*pow(physcon.pi,2)
//                                                           *physcon.epsilon_0*
//                                           physcon.c*physcon.m_e)
//     *kappa*BField*pow(((3*physcon.e)/(2*physcon.pi*nu*pow(physcon.m_e,3)*
//                          pow(physcon.c,4))),p/2.)*b_p*pow(nu,-temp2);
  
 //  std::cout<<"abs"<<"\t"<<ChiNu<<"\t"<<pow(BField,temp1)
//            <<"\t"<<pow(nu,temp2)<<"\t"<<exp(-ChiNu)<<std::endl;
   
  //units m^-1
  return ChiNu;
  
}

//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
