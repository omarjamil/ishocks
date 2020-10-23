//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....//
// statistical analysis function such as Lomb-Scargle                  //
// Omar Jamil 2008                                                     //
//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....//

#ifndef STATANALYSIS_HH
#define STATANALYSIS_HH

#include <vector>
#include <cmath>
#include <string>

//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....

class StatAnalysis
{
  
public:
  StatAnalysis(); 	//the default constructor
  //statAnalysis(const statAnalysis &c);
  ~StatAnalysis();
  //protected:

  void aveVar(std::vector<double> &, double &, double &);
  
  void periodogram(std::vector<double> &, 
                   std::vector<double> &, 
                   const double, const double, 
                   std::vector<double> &, std::vector<double> &,
                   double &, int &, double &);

  void getData(std::string &, std::vector<double> &, 
               std::vector<double> &);
  
  void writeData(std::string &, std::vector<double> &,
                 std::vector<double> &);

  void writeData(std::string &,
                 std::vector<double> &);
  
  void fasPer(std::vector<double> &, 
              std::vector<double> &,
              const double, const double,
              std::vector<double> &, 
              std::vector<double> &,
              double &, int &, double &
              );
  
  void spread(const double , std::vector<double> &, 
              const double , const int );
  
  void realFT(std::vector<double> &, const int);
  
  void four1(std::vector<double> &, const int);

  void fourierOnly(std::vector<double> &, const int);
  
   
    
  //private:
  
};

//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
class Parameters
{
  
public:
  Parameters(std::string); 	//the default constructor
  //Parameters(const Parameters &c);
  ~Parameters();
  
  void loadParametersFile(std::string) throw (const std::exception &);

   
  std::string inputFile()
  {
    return inputFile_;
  }

  std::string outputFile()
  {
    return outputFile_;
  }
  
  int analysisChoice()
  {
    return analysisChoice_;
  }
  
  double overSampling()
  {
    return overSampling_;
  }
  
  double freqFactor()
  {
    return freqFactor_;
  }
  
  
     
private:
  std::string inputFile_;
  std::string outputFile_;
  int analysisChoice_;
  double overSampling_;
  double freqFactor_;
  
  
};
//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....



#endif // STATANALYSIS_HH
