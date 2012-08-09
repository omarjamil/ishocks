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
#include <iterator>
#include <fstream>

#include "shells.hh"
#include "physcon.hh"
#include "parameters.hh"
#include "paraccess.hh"

//a small macro for square function
#define sqr(a) ((a)*(a))

//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....

//constructor with member initialize
//present in .cc in case more complicated initialization req.
Shells::Shells(double &i, double &m, double &g, double &e, double &w)
  : ti(i), mass(m), gamma(g), ienergy(e), width(w), betaExpansion(0.0),
    x_l(0.0), x_u(0.0), r_l(0.0), r_u(0.0), vol(0.0), area(0.0),
    egammaMin(exteGammaMin), egammaMax(exteGammaMax), mergeCount(0),
    shell(true), pNorm(0.0), BEDens(0.0), magPressure(0.0), shellCentre(0.0),
    tPrev(i), tEnergize(0.0), fracEInt(0.0), intEavail(0.0), shocked(0)
{
  try
    {
      iNuGrid = new std::vector<double>;
      tau = new std::vector<double>;
      gammaGrid = new std::vector<double>;
      gammaRange();

    }
  catch(std::bad_alloc)
    {
      std::cerr<<"Memory exhauseted in constructor;shell.cc\n";
    }

  //Various shell initialization upon creation
  shellThermEner = extShThermalEne * (e / g);

}

//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
//Compiler Copy constructor in use
//Shells::Shells(const Shells&): tc(c),{}
//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....


//destructor
Shells::~Shells()
{
  delete iNuGrid;
  delete tau;
  delete gammaGrid;

}

//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
void Shells::setLowerAndUpperX(double &LowerX, double &UpperX)
{
  setInnerRadius(LowerX);
  setOuterRadius(UpperX);

}
//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
void Shells::setMergeCount(int &totalCount)
{
  mergeCount = totalCount + 1;
}
//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
void Shells::shellInitialization()
{
  double wid = width;
  double theta = (extJetOpenAngle/180.0)*physcon.pi; //jet opening angle

  x_l = 1.0;
  //x_u = 0.5*(wid);
  x_u = wid;

  r_l = 2.* x_l * tan(theta/2.);
  r_u = 2.* x_u * tan(theta/2.);

  vol = (1.0/3.0) * physcon.pi * (0.5*wid) *
    (sqr(r_l) + (r_l * r_u) + sqr(r_u));


  area = physcon.pi*(r_u*sqrt(sqr(r_u) + sqr(x_u)) -
                     r_l*sqrt(sqr(r_l) + sqr(x_l)) );



}

//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
void Shells::setDimensions(double &t_ela)
{

  //shell width (initial)
  double betaExp = getExpansionBeta();
  //double wid = width + (betaExp * physcon.c * (t_ela - ti));
  //std::cout<<width<<std::endl;

  width += (betaExp * physcon.c * (t_ela - tPrev));
  //width = wid;

  //shell velocity
  double beta = sqrt(1. - (1./sqr(gamma)));

  shellCentre = beta * physcon.c * (t_ela - ti);

  //setting x_lower and x_upper
  // x_l = beta * physcon.c * (t_ela - ti) - (betaExp * physcon.c *
//                                            (t_ela - ti));

//   x_u = beta * physcon.c * (t_ela - ti) + (betaExp * physcon.c *
//                                            (t_ela - ti));

  x_l = shellCentre - 0.5*(width);
  if(x_l < 0.0)
    {
      x_l = 1.0;
    }
  x_u = shellCentre + 0.5*(width);


  //std::cout<<x_l<<"\t"<<x_u<<"\t"<<"x_l and x_u"<<std::endl;

  double theta = (extJetOpenAngle/180.0)*physcon.pi; //jet opening angle
  //Shell radius at lower x value
  r_l = 2.* x_l * tan(theta/2.);
  r_u = 2.* x_u * tan(theta/2.);

  double shWid = width;

  //conical frustum volume
  vol = (1.0/3.0) * physcon.pi * shWid *
    (sqr(r_l) + (r_l * r_u) + sqr(r_u));
  //std::cout<<t_ela<<"\t"<<ti<<"\t"<<vol<<"\t"<<shWid<<"\t"<<"in setDim"<<std::endl;
  //std::cout<<ti<<"\t"<<betaExp<<"\t"<<t_ela-tPrev<<"\t"<<vol<<"\t"<<"in setDim"<<std::endl;


    //conical frustum surface area
 //  area = physcon.pi*(r_u*sqrt(sqr(r_u) + sqr(x_u)) -
//                      r_l*sqrt(sqr(r_l) + sqr(x_l)) );

  area = physcon.pi*(r_u+r_l)*sqrt(sqr(r_u-r_l)+sqr(shWid));

  tPrev = t_ela;

  shWid = 0.0;

}
//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
void Shells::setInitialMergedVA()
{
  double theta = (extJetOpenAngle/180.0)*physcon.pi;
  double wid = width;

  r_l = 2.* x_l * tan(theta/2.);
  r_u = 2.* x_u * tan(theta/2.);

  vol = (1.0/3.0) * physcon.pi * wid *
    (sqr(r_l) + (r_l * r_u) + sqr(r_u));



  //and area
 //  area = physcon.pi*(r_u*sqrt(sqr(r_u) + sqr(x_u)) -
//                      r_l*sqrt(sqr(r_l) + sqr(x_l)));

  area = physcon.pi*(r_u+r_l)*sqrt(sqr(r_u-r_l)+sqr(wid));

  //std::cout<<r_l<<"\t"<<r_u<<"\t"<<vol<<"\t"<<area<<std::endl;


}

//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....

void Shells::setExpansionBeta()
{
  //Spada et al '01 equations 7 and 8
  //fraction of internal energy given to thermal energy
  double shellMass = getShellMass();
  double shellGamma = getShellGamma();
  //double iEnergy = getInternalEnergy() / shellGamma;
  double thermalEnergy = getThermalEnergy();
  double soundVelocity = sqrt((1./3.)*(thermalEnergy/shellMass));
  double betaS = soundVelocity/physcon.c;
  double shellbeta = sqrt(1. - (1. / sqr(shellGamma)));

  betaExpansion = ((2.*betaS)/sqr(gamma))*
     (1./(1.-pow((shellbeta*betaS),2)));


  //betaExpansion = 0.0;


}

//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
void Shells::powerlawNorm()
{
  double shellVolume = vol;//this->getShellVolume();
  double /*kappa = 0.0,*/ kappa2=0.0;
  double plindex = extPowLawInd;
  //N(E) in GeV for calculating radiation - Longair
  double iEnergy = getInternalEnergy();
  //6.241590974e18 eV = 1J
  //the following in units of GeV; divided by 511e-6 eV as
  //the calculation is done in units of mc2
  double EintK = (((iEnergy*extShIntKine)/shellVolume)*6.241590974e9)/(511e-6);
  double EintK2 = ((iEnergy*extShIntKine)/shellVolume)/(physcon.mc2());
  double gammaMax = getGammaMax();
  double gammaMin = getGammaMin();
  // double newIntEn = iEnergy - (iEnergy*extShIntKine);
//   setInternalEnergy(newIntEn);

  if(EintK == 0.0)
    {
      //kappa = 0.0;
      kappa2 = 0.0;
    }
  else if(gammaMax == 1.)
    {
      //kappa = 0.0;
      kappa2 = 0.0;
    }
  else
    {

      if(plindex == 2.0)
        {
	  //  kappa = pow(511e-6, (plindex-1.)) *
          //  EintK/((log(gammaMax) - log(gammaMin)) + ((1./(gammaMax))
          //                                            - (1./(gammaMin))));

	  kappa2 = pow(physcon.mc2(), (plindex-1.)) *
            EintK/((log(gammaMax) - log(gammaMin)) + ((1./(gammaMax))
                                                      - (1./(gammaMin))));

	  

        }
      else
        {

	  // kappa = pow(511e-6, (plindex-1.)) *
          //  EintK/(((1./(2.-plindex))*(pow(gammaMax,(2.-plindex)) -
          //                             pow(gammaMin,(2.-plindex)))) -
          //         ((1./(1.-plindex))*(pow(gammaMax,(1.-plindex)) -
          //                             pow(gammaMin,(1.-plindex)))));

          //use with SI units calc. of synchrotron
          kappa2 = pow(physcon.mc2(), (plindex-1.)) *
            EintK2/(((1./(2.-plindex))*(pow(gammaMax,(2.-plindex)) -
                                       pow(gammaMin,(2.-plindex)))) -
                   ((1./(1.-plindex))*(pow(gammaMax,(1.-plindex)) -
                                       pow(gammaMin,(1.-plindex)))));

        }
    }

  pNorm = kappa2;
  powerLaw();


}
//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
void Shells::adiabaticLosses(double &tNow)
{
  bool shell = this->getShellId();

  if(!shell)

    {
      double shVolume = getShellVolume(); //Shell volume before update

      double plind = extPowLawInd;

      //powerlaw normalization before update
      double pNormBefore = getPLawNorm();

      double shVolumeNow, shWidth, shbeta, xLower, xUpper, rLower,
        rUpper, shTheta, pNormAfter, gammaMaxBefore, gammaMaxAfter,
        magPressureAfter, magPressureBefore, intEneBefore, intEneAfter,
        /*thermalEneBefore,*/ thermalEneAfter;

      double betaExp = this->getExpansionBeta();
      double shellGamma = this->getShellGamma();

      gammaMaxBefore = this->getGammaMax();
      intEneBefore = this->getInternalEnergy();


      shWidth = width + (betaExp * physcon.c * (tNow - tPrev));
      shbeta = sqrt(1. - (1./sqr(gamma)));
      // xLower = shbeta * physcon.c * (tNow - ti) -
//         (betaExp * physcon.c * (tNow - ti));
//       xUpper = shbeta * physcon.c * (tNow - ti) + shWidth;

      xLower = (shbeta * physcon.c * (tNow - ti)) - (0.5 * shWidth);
      if(xLower < 0.0)
        {
          xLower = 1.0;
        }
      xUpper = (shbeta * physcon.c * (tNow - ti)) + (0.5 * shWidth);

      shTheta = (extJetOpenAngle/180.0)*physcon.pi; //jet opening angle
      rLower = 2.* xLower * tan(shTheta/2.);
      rUpper = 2.* xUpper * tan(shTheta/2.);

      //what shell volume will be after update
      shVolumeNow = (1.0/3.0) * physcon.pi * shWidth *
        (sqr(rLower) + (rLower * rUpper) + sqr(rUpper));
      //std::cout<<tNow<<"\t"<<ti<<"\t"<<width<<"\t"<<shVolume<<"\t"<<shVolumeNow<<std::endl;

      //std::cout<<ti<<"\t"<<betaExp<<"\t"<<tNow-tPrev<<"\t"<<shVolumeNow<<"\t"<<"in adia"<<std::endl;

      //std::cout<<ti<<"\t"<<shVolumeNow<<"\t"<<shVolume<<"\t"<<tNow-ti<<"\t"<<shWidth<<" in adia"<<std::endl;
      double temp = (-plind-2.)/3.;

      //the new powerlaw normalization
      pNormAfter = pNormBefore * pow((shVolumeNow/shVolume), temp);
      this->setPLawNorm(pNormAfter);
      powerLaw();

      double temp2 = -1./3.;
      //New gammaMax after adiabatic losses;
      gammaMaxAfter = gammaMaxBefore * pow((shVolumeNow/shVolume), temp2);
      if(gammaMaxAfter <= 1.0)
        {
          gammaMaxAfter = 1.0;
        }
      //std::cout<<gammaMaxBefore<<"\t"<<gammaMaxAfter<<"\t"<<pow((shVolumeNow/shVolume), temp2)<<std::endl;

      setGammaMax(gammaMaxAfter);

      //thermalEneBefore = getThermalEnergy();
      intEneAfter = intEneBefore * pow((shVolumeNow/shVolume), temp2);
      this->setInternalEnergy(intEneAfter);
      //thermalEneAfter = thermalEneBefore * pow((shVolumeNow/shVolume), temp2);
      thermalEneAfter = extShThermalEne * (intEneAfter / shellGamma);
      setThermalEnergy(thermalEneAfter);
      //this->setExpansionBeta();

      pNormBefore = 0.0;
      double temp3 = -4./3.;
      magPressureBefore = getMagPressure();
      magPressureAfter = pow((shVolumeNow/shVolume), temp3)
        * magPressureBefore;
      //double u = getBEneDens();
      this->setMagParams(magPressureAfter);
      //double v = getBEneDens();

      //std::cout<<ti<<"\t"<<magPressureBefore<<"\t"<<magPressureAfter<<"\t"<<"Mag"<<std::endl;


    }
  else
    {
      pNorm = 0.0;
      powerLaw();
    }


}
//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
double Shells::dopplerFactor()
{

  double viewAngle = extJetViewAngle/180. * physcon.pi;
  double shellGamma = gamma;
  double shellBeta = sqrt(1. - (1. / sqr(gamma)));

  //approaching jet
  double dF = 1./(shellGamma*(1.-shellBeta*cos(viewAngle)));

  //receding jet
  //double dF = 1./(shellGamma*(1.+shellBeta*cos(viewAngle)));

 return dF;


}
//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
void Shells::initialBEnergyDensity()
{
  //calculated in the shell rest frame.

  double shellVolume = getShellVolume();
  double intEne = getInternalEnergy();

  BEDens =  sqrt(2.*(((intEne*extShIntMag)/shellVolume)*physcon.mu_0));
  magPressure = BEDens;

  // double newIntEn = intEne - (intEne*extShIntMag);
//   setInternalEnergy(newIntEn);
}
//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
void Shells::setMagParams(double &newMagPressure)
{
  if(newMagPressure <= 0.0)
    {
      magPressure = 0.0;
    }
  else
    {
      magPressure = newMagPressure;
    }

  BEDens = newMagPressure;

}
//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
void Shells::inVacuumExpansion()
{
  bool shell = this->getShellId();
  if(!shell)
    {
      this->powerlawNorm(); //powerlaw norm for bigger shell
      this->initialBEnergyDensity(); //Magnetic energy density.
      //thermal energy remains constant.(probably should be kept zero in par file).
      double thermEne = getThermalEnergy();
      this->setThermalEnergy(thermEne);

    }

}
//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
void Shells::slowEnergization(double &tNow)
{
  bool shell = this->getShellId();

  if(!shell)
    {
      double tEnerg = this->getTEnergize();
      double tRatio = tNow/tEnerg;

      double prevEIntFrac = this->getIntEfrac();
      double intEneShell = this->getInternalEnergy();
      double intEneAvail = this->getIntEavail();
      double intEfracNow = intEneAvail*tRatio;
      double EShell = intEneShell + (intEfracNow - prevEIntFrac);

      if(tRatio <= 1)
        {
          //std::cout<<tNow<<"\t"<<tEnerg<<"\t"<<tEnerg<<"\t"<<EShell<<std::endl;
          this->setInternalEnergy(EShell);
          this->setIntEfrac(intEfracNow);
          this->powerlawNorm(); //re-calc powerlaw norm with increased
                                //internal energy
          this->initialBEnergyDensity(); //re-calc B energy density

        }
      //if the shell has some internal energy and past
      //tFrac = 1.0 keep the status quo i.e. only the evolution
      else if(tRatio > 1.0 && intEneShell > 0.0)
        {

          this->setInternalEnergy(intEneShell);


        }
      //once past tFrac = 1.0 give all the internal energy
      //to the shell
      else if(tRatio > 1.0 && intEneShell == 0.0)
        {
          this->setInternalEnergy(intEneAvail);
          this->powerlawNorm(); //re-calc powerlaw norm with increased
                                //internal energy
          this->initialBEnergyDensity(); //re-calc B energy density
          //std::cout<<tNow<<"Int En == 0.0"<<std::endl;
        }
    }

}
//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
void Shells::shockZone(double &tNow)
{
  double shellLoc = this->getLocation();
  double injT = this->getTimeOfInjection();
  double relMassFrac = extRelMassFrac;
  double shellG = this->getShellGamma();
  double shellM = this->getShellMass();
  double intEneSet = 0.0;
  bool shShock = this->shellShocked();
  double shLoc = extShockLoc;
  double shocLoc = shLoc * physcon.c;

  //the following is done because if shocLoc zero the shell has internal
  //energy at location=0.0, but very small volume (initial size).
  //The next time step the shell gets the proper volume which is
  //much larger than the initial volume. The relative change causes the
  //shell to lose all it's internal energy very quickly. By setting shocLoc at only
  //a few meters, the shell does not get internal energy until the following
  //time step where it has achieved it's normal volume first. The order of
  //things done in Container::updateShellParams causes this discrepancy to
  //arise. Shock Location of a few metres ~ 0.0!
  if(shocLoc == 0.0)
    {
      shocLoc = 3.e-8;
    }

  if(shellLoc >= shocLoc)
    {
      //std::cout<<tNow<<"\t"<<injT<<"\t"<<shellLoc<<"\t"<<shocLoc<<"\t"<<intE<<std::endl;

      if(!shShock)
        {
         intEneSet = relMassFrac * shellG * shellM * (physcon.c * physcon.c);
         //intEneSet = relMassFrac * shellM * (physcon.c * physcon.c);
         this->setTEnergize(injT);
         this->setInternalEnergy(intEneSet);
         this->setIdFalse();
         this->powerlawNorm();
         this->setExpansionBeta();
         this->initialBEnergyDensity();
         this->shockShell();
        }
    }

}
//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
void Shells::radiativeLosses(double &tNow)
{
  bool shell = this->getShellId();
  double TCS = 6.6524586e-29;
  double eMass  = 9.109382e-31;
  double cSpeed = 3.e8;
  double mu_0 = 1.257e-6;

  if(!shell)
    {
      //double gammaMaxBefore = this->getGammaMax();
      double magP = this->getMagPressure();
      double intEneBefore = this->getInternalEnergy();

      double cons = 1.5*eMass*eMass*cSpeed*cSpeed*cSpeed*mu_0/TCS;

      double eneNow = intEneBefore/(1.+ (intEneBefore * cons
                                         * magP * magP * (tNow-tPrev)));


      //std::cout<<cons<<"\t"<<magP<<std::endl;

      // std::cout<<eneNow<<"\t"<<intEneBefore<<std::endl;
      // std::cout<<intEneBefore-eneNow<<std::endl;

      this->setInternalEnergy(eneNow);
      this->powerlawNorm();
      this->setExpansionBeta();
      this->initialBEnergyDensity();


    }
}
//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
void Shells::gammaRange()
{
  egamma = new std::vector<double>;
  double min = egammaMin;
  double max = egammaMax;

  int gammaPoints = extgammaPoints;

  double abscissa = min;
  double r = pow(max / min, 1.0 / gammaPoints);
  double r2 = sqrt(r);
  r -= 1.0;

  //writing the frequency range file
  //useful when sampling many frequencies

  //std::ofstream gammaFile("gammaRange.dat");
  std::cout<<"Log binning for electrons(Lorentz factor)"
	   <<std::endl;
  for (int i = 1; i <= gammaPoints; ++i)
    {
      egamma->push_back(r2 * abscissa); // Calc log-space mid-point for grid
      abscissa += (r * abscissa);  // Add grid width to abscissa to find
      //std::cout<<r2<<"\t"<<(abscissa)<<std::endl;
      //gammaFile<<r2 * abscissa<<std::endl;

    }

}

//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
void Shells::powerLaw()
{
  //double n0 = ne*(p-1.0);
  double n0 = getPLawNorm();
  double gammaMin = getGammaMin();
  double gammaMax = getGammaMax();

  for (std::vector<double>::iterator n = gammaGrid->begin(),
	 gamma = egamma->begin(); gamma != egamma->end(); ++n, ++gamma)
    {
      // quicker than starting at begin()+index(gamma_min)?
      if (*gamma >= gammaMin)
	{
	  if (*gamma <= gammaMax)
	    {
	      *n += n0 * pow(*gamma, -extPowLawInd);
	    }
	  else
	    {
	      *n = 0.0;
	    }
	}
    }
}

//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
