# ifndef __ADAMS_METHOD_SOLVERS_H__
# define __ADAMS_METHOD_SOLVERS_H__ 

# include "../../rhsODEproblem.H"
# include "../MultiStep.H"


namespace mg { 
                namespace numeric {
                                    namespace odesystem {


/*-------------------------------------------------------------------------------
 *    
 *    Base class for Adams Methods : 
 *    
 *    (*) Adams Bashforth (2nd,3th order) explicit method
 *    (*) Adams Moulton (2nd , 3th order) Predictor Corrector (implicit method)
 *    
 *
 *    @Marco Ghiani Dec 2017, Glasgow UK
 *
 ------------------------------------------------------------------------------*/


template <typename Type > 
class AdamsMethods            :      
                                  public MultiStep<Type>
{
    
   public: 
      
      AdamsMethods(const rhsODEProblem<Type>& that ) noexcept : 
                                                              MultiStep<Type>{that} 
                          {}                                
      
      virtual ~AdamsMethods() = default ;
      
      using OdeSolver<Type>::rhs;
      
      virtual void solve(const std::string filename) override = 0;
      virtual void solve() noexcept  override                 = 0;


   protected:
     
      std::valarray<Type> um1;
      std::valarray<Type> um2;
      std::valarray<Type> um3;
      std::valarray<Type> um4;
      std::valarray<Type> up1;
      std::valarray<Type> u;

        
      
      std::valarray<Type> K1 ;
      std::valarray<Type> K2 ;
      std::valarray<Type> K3 ;
      std::valarray<Type> K4 ;
      std::valarray<Type> K5 ;
      std::valarray<Type> K6 ;
      
      std::valarray<Type> uPred ; // u predictor
      std::valarray<Type> uCorr    ; // u corrector
      std::valarray<Type> uCorrOld ; // u corrector

      std::valarray<Type> fPred ;    // f(time+1) predictor
      std::valarray<Type> fCorr ;    // f(time+1) corr
      std::valarray<Type> fCorrOld ;    // f(time+1) corr
      
      std::valarray<Type> error ;    // iteration check

};



  }//ode 
 }//numeric
}//mg
# endif 
