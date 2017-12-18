# ifndef __ABSTRACT_ODE_SOLVER_H__
# define __ABSTRACT_ODE_SOLVER_H__

# include <exception>
# include <string>
# include <vector>
# include <valarray>
namespace mg {
               namespace numeric { 
                                    namespace odesystem {

template<typename Type>
class AbstractODESolver
{
    public:
        
        virtual void setStepSize() = 0 ;
        virtual void setInitialTime() = 0 ;
        virtual void setFinalTime () = 0 ;
        virtual void setInitialValue() =0 ;

        virtual Type getStepSize()    const = 0 ;
        virtual Type getInitialTime() const = 0 ;
        virtual Type getFinalTime()   const = 0;
        virtual std::valarray<Type> getInitialValue()const = 0 ;
        
        virtual void solve(const std::string filename)  = 0;
        virtual void solve() noexcept             = 0;
        virtual ~AbstractODESolver() = default;



   protected:
  /*   
   
   */   
};

  }//ode
 }//numeric
}//mg
# endif
