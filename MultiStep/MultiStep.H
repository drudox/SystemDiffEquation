# ifndef __MULTI_STEP_SOLVERS_H__
# define __MULTI_STEP_SOLVERS_H__ 

# include "../rhsODEproblem.H"
# include "../OdeSolver.H"


namespace mg { 
                namespace numeric {
                                    namespace odesystem {


/*-------------------------------------------------------------------------------
 *    
 *    Base class for Multi-Step solver of differential IVP
 *
 *
 *    @Marco Ghiani Dec 2017, Glasgow UK
 *
 ------------------------------------------------------------------------------*/


template <typename Type > 
class MultiStep            :      
                                  public OdeSolver<Type>
{
    
   public: 
      
      MultiStep(const rhsODEProblem<Type>& that ) noexcept : 
                                                              OdeSolver<Type>{that} 
                          {}                                
      
      virtual ~MultiStep() = default ;
      
      using OdeSolver<Type>::rhs;
      
      virtual void solve(const std::string filename) override = 0;
      virtual void solve() noexcept  override                 = 0;


   protected:
     
      std::valarray<Type> u_m1;
      std::valarray<Type> u_p1;
      std::valarray<Type> u_  ;
      
      std::valarray<Type> k1;
      std::valarray<Type> k2;

};



  }//ode 
 }//numeric
}//mg
# endif 
