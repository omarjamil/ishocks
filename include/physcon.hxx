//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo....//
//                                                                    //
// Internal shock model for relativistic astrophysical jets           //
//                                                                    //
// Author: Omar Jamil                                                 //
// Contact: o.jamil@phys.soton.ac.uk                                  //
//          University of Southampton                                 //
//                                                                    //
//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo....//

#ifndef PHYSCON_HH
#define PHYSCON_HH

// Precision is at least 7 sig fig

class PhysicalConstants {
public:
  PhysicalConstants() : pi(3.14159265359), 
                        c(2.99792458e8 /*m s^-1*/), 
                        e(1.60217733e-19), 
                        m_e(9.1093897e-31), 
                        sigma_t(6.652462e-29), 
                        h(6.6260755e-34 /*Js*/), 
                        alpha_f(7.29735308e-3), 
                        epsilon_0(8.8541878e-12), 
                        k(1.380658e-23 /*J K^-1*/), 
                        pc(3.085678e16), 
                        sigma(5.67040040e-8/*W m^-2 K^-4*/), 
                        h_keV(4.135667e-18 /*keV s*/),
                        mu_0(12.56637061e-7 /*N A^-2*/)
  { }
  double mc2() const { return m_e * c * c; }
  const double pi;
  const double c;
  const double e;
  const double m_e;
  const double sigma_t;
  const double h;
  const double alpha_f;
  const double epsilon_0;
  const double k;
  const double pc;
  const double sigma;
  const double h_keV;
  const double mu_0;
  
};

const PhysicalConstants physcon;

#endif //PHYSCON_HH
