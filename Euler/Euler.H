# ifndef __EULER_SOLVERS_INTERFACE_H__
# define __EULER_SOLVERS_INTERFACE_H__

# include "../OdeSolver.H"
# include "../rhsODEproblem.H"

namespace mg {
                namespace numeric {
                                     namespace odesystem {


/*-----------------------------------------------------------------------------
 *       
 *    @class interface for Eulers method (implicit and explicit)  
 *           solvers of ODE (RHS problem) IVP 
 *
 *    @autor Marco Ghiani 
 *    @date  Dec 2017, Glasgow UK
 *
 -----------------------------------------------------------------------------*/




template<typename Type = double>
class Euler :
                           public OdeSolver<Type> 
{
      
    public:  
      Euler(const rhsODEProblem<Type> & that) noexcept :
                                                                        OdeSolver<Type>{that} 
                  {}
      
      virtual ~Euler() = default ;

      using OdeSolver<Type>::rhs;

      void solve(const std::string filename) override = 0  ;
      void solve() noexcept override  = 0                  ;
//
//
      unsigned short order() {return 1; } // return the order of the solvers 
};

 
  }//ode
 }//numeric
}//mg
# endif
