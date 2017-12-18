# ifndef __ODE_SOLVER_H__
# define __ODE_SOLVER_H__

# include "abstractOdeSolver.H"
# include <cmath>
# include <iostream>
# include <cstdlib>
# include "rhsODEproblem.H"
# include <vector>
# include <valarray>
# include <string>
//# include "rhsOdeProblem.H"
//# include "RHS_ODE.H"

namespace mg { 
                namespace numeric {
                                     namespace odesystem {

template< typename T>
class rhsODEProblem;


/*------------------------------------------------------------------------*
 *    
 *    Abstract OdeSolver Class - long double double or single precision  
 *    
 *    @author Marco Ghiani 
 *    @date   Dec 2017 Glasgow
 *
 -------------------------------------------------------------------------*/


template <typename Type>
class OdeSolver  : 
                    public AbstractODESolver<Type> 
{
  
//--
//
  public:

    OdeSolver(const rhsODEProblem<Type>& that) noexcept : rhs{that} 
    {
     
     setStepSize    () ; 
     setInitialTime () ;
     setFinalTime   () ;
     setInitialValue() ;
     setSize() ;
    }

    virtual ~OdeSolver() = 0;

    rhsODEProblem<Type> rhs ;
    

     virtual void setStepSize    () override { stepSize     = rhs.dt(); }
     virtual void setInitialTime () override { initialTime  = rhs.t0(); }
     virtual void setFinalTime   () override { finalTime    = rhs.tf(); }
     virtual void setInitialValue() override 
     {
            initialValue = rhs.u0(); 
     }

     virtual void solve(const std::string filename) = 0;
     virtual void solve() noexcept            = 0;
             auto setProblem(const rhsODEProblem<Type>& that ) { this->rhs = that ; }  
     
     virtual Type getStepSize()    const override { return stepSize     ;}
     virtual Type getInitialTime() const override { return initialTime  ;}
     virtual Type getFinalTime()   const override { return finalTime    ;}
     virtual std::valarray<Type> getInitialValue() const override { return initialValue ;}
      
     virtual Type dt() const noexcept { return stepSize     ;}
     virtual Type t0() const noexcept { return initialTime  ;}
     virtual Type tf() const noexcept { return finalTime    ;}
     virtual std::valarray<Type> u0() const noexcept { return initialValue ;}
     
     virtual void setSize() noexcept ;
     
     virtual void setTimeStep(const Type& h_ ) noexcept {   // used from the adaptative Time step solver
          stepSize = h_ ;  
     }

     protected:
      
      Type  stepSize;
      Type  initialTime;
      Type  finalTime;
      std::valarray<Type> initialValue;
      
      const int Ns = (rhs.tf() - rhs.t0())/rhs.dt() ;
      
      constexpr static Type toll = 1e-12 ;
      
      Type                t ;
      std::valarray<Type> u ;
      
};

template<typename Type>
OdeSolver<Type>::~OdeSolver() = default ;

template<typename Type>
void OdeSolver<Type>::setSize() noexcept
{
  u.resize(initialValue.size()) ;
}


  }//ode
 }//numeric
}//mg
# endif
