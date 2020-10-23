//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....//
// statistical analysis function such as Lomb-Scargle                  //
// Omar Jamil 2008                                                     //
//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....//

#include <vector>
#include <iostream>

#include "statAnalysis.hh"


int main()
{
  Parameters params("parameters.par");
  
  std::vector<double> x;
  std::vector<double> y;
  
  std::string file = params.inputFile();
  std::string outFile = params.outputFile();
  
  double ofac = params.overSampling();
  double hifac = params.freqFactor();
    
  int choice = params.analysisChoice();

  double freqMax;
  double nout = 100;
  int jmax = 0;
  double prob = 0.0;
  
  StatAnalysis stats;

  stats.getData(file,x,y);
  int xsize = x.size();
  //std::cout<<xsize<<std::endl;
  
  double fc = x.size()/(2*(x.back() - x.front()));
  
  if (choice != 3)
    {
      std::cout<<"The maximum frequency is: "<<hifac * fc<<std::endl;
    }
  
  double np_double = 0.5*xsize*ofac*int(hifac+1.);
  int np = int(np_double);
   

  std::vector<double> px(np);
  std::vector<double> py(np);

  std::vector<double> wk1(np);
  std::vector<double> wk2(np);

  if(choice == 2)
    {
      stats.periodogram(x,y,ofac,hifac,px,py,nout,jmax,prob);
      stats.writeData(outFile, px, py);
      
      std::cout<<"Number of frequencies: "<<nout<<std::endl;
      std::cout<<"Maximum at freq: "<<px[jmax]<<std::endl;
      std::cout<<"Null Probability: "<<prob<<std::endl;
    }
  else if(choice == 1)
    {
      
      stats.fasPer(x,y,ofac,hifac,wk1,wk2,nout,jmax,prob);
      stats.writeData(outFile,wk1,wk2);
    
      std::cout<<"Number of frequencies: "<<nout<<std::endl;
      std::cout<<"Maximum at freq: "<<wk1[jmax]<<std::endl;
      std::cout<<"Null Probability: "<<prob<<std::endl;
    }
  else if(choice == 3)
    {
      stats.fourierOnly(y,1);
      stats.writeData(outFile,y);
      std::cout<<"Read y data and wrote 2 column file"<<"\n"; 
      std::cout<<"Column 1 = real; Column 2 = imaginary"
               <<std::endl;
      
    }
  
      
  return 1;
  
  
}

