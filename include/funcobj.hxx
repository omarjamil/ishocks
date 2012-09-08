//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo....//
//                                                                    //
// Internal shock model for relativistic astrophysical jets           //
//                                                                    //
// Author: Omar Jamil                                                 //
// Contact: o.jamil@phys.soton.ac.uk                                  //
//          University of Southampton                                 //
//                                                                    //
//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo....//

#ifndef FUNCOBJ_HH
#define FUNCOBJ_HH


//collection of function objects etc.
// The following function objects are used as predicates for
// various finding/sorting/inserting algorithms provided by STL

#include <functional>
#include <cmath>


#include "shells.hxx"
#include "mergers.hxx"

//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....

class FindShellWithShell
{
  Shells *shell;
  
  
public:
  FindShellWithShell(Shells *s) : shell(s){}
  

  bool operator()(Shells *data)
  {
    return (shell->getTimeOfInjection() == data->getTimeOfInjection());
    
  }
  
  
};

//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....

class FindShellWithDouble
{
  double t;
       
public:
  FindShellWithDouble(double &s) : t(s){}
   

  bool operator()(Shells *data)
  {
    return (data->getTimeOfInjection() == t);
     
  }
  
  
};
//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....


class FindLowerBound
{
  //double t;
       
public:
  FindLowerBound(){}
   

  bool operator()(Shells *data, const double &t)
  {
    return (data->getTimeOfInjection() < t);
     
  }
  
  
};

//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....

class FindUpperBound
{
  //double t;
       
public:
  FindUpperBound(){}
   

  bool operator()(const double &t, Shells *data)
  {
    return (data->getTimeOfInjection() > t);
     
  }
  
  
};

//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....


class XFindLowerBound
{
  //double t;
       
public:
  XFindLowerBound(){}
   

  bool operator()(Shells *data, const double &t)
  {
    return (data->getInnerRadius() > t);
     
  }
  
  
};

//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
// Function object structure that inherits from STL binary
// functions. Return operator () defined to return a result.
// Only a bool required as this is used with STL 'stable_sort' as 
// 'comp' function. Structure used as only a simple thing with only
// public members required, no initialization. No arguments being
// passed. Class can be used with public decl only.


//The following cannot be used with list as sort uses random
//access iterators which are not available for std::list<>

class Sortfunc : public std::binary_function<Shells*, Shells*, bool>
{
  
public:  
  
  bool operator()(Shells* lhs, Shells* rhs)
    {
      // use '<' to sort ascending 
      // use '>' to sort descending
      return lhs->getTimeOfInjection() < rhs->getTimeOfInjection();
    }
};

//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
class deleteAll
{
 public:
  deleteAll(){}
   
  bool operator()(Mergers *theElement)
  {
    delete theElement; 
    return true;
  }
};

//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
template< typename T >
struct delete_ptr : public std::unary_function<bool,T>
{
  bool operator()(T*pT) 
    const 
  { 
    delete pT; 
    return true; 
  }
};
//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....



#endif // FUNCOBJ_HH
