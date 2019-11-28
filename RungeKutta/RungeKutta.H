# ifndef __RUNGEKUTTA_H__
# define __RUNGEKUTTA_H__ 

# include "../rhsODEproblem.H"
# include "../OdeSolver.H"


namespace mg { 
                namespace numeric {
                                    namespace odesystem {


/*-------------------------------------------------------------------------------
 *    
 *
 *    @class Base class for all the Runge-Kutta methods 
 *    dy/dt = f(y,t)
 *
 *    @author Marco Ghiani October 2017, Glasgow UK
 *
 ------------------------------------------------------------------------------*/


template <typename Type = double> 
class RungeKutta :      
                        public OdeSolver<Type>
{
    
   public: 
      
      RungeKutta(const rhsODEProblem<Type>& that ) noexcept : 
                                                          OdeSolver<Type>{that} 
                          {}                                
      
      virtual ~RungeKutta() = default ;
      
      using OdeSolver<Type>::rhs;
      
    //  virtual void solve(const std::string& ) override = 0 ;
    //  virtual void solve() noexcept override           = 0 ;
      virtual void solve(const std::string filename) override = 0;
      virtual void solve() noexcept  override                 = 0;


   protected:
     
      std::valarray<Type> k1;
      std::valarray<Type> k2;
      std::valarray<Type> k3;
      std::valarray<Type> k4;
      std::valarray<Type> k5;
      std::valarray<Type> k6; 
      std::valarray<Type> k7;
      std::valarray<Type> k8; 
      
      std::valarray<Type> up; 
      std::valarray<Type> uc; 
      std::valarray<Type> u;
      std::valarray<Type> u_5;
      std::valarray<Type> u_4;
      
      std::valarray<Type> s ;
      std::valarray<Type> truncError ;
      
      const static Type dt0 ;

      constexpr static Type stepToll = 1.0e-6 ;
      
};



  }//ode 
 }//numeric
}//mg
# endif 
