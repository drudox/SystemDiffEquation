# ifndef __FORWARD_EULER_SOLVER_H__
# define __FORWARD_EULER_SOLVER_H__

# include "Euler.H"
# include "../rhsODEproblem.H"

namespace mg {
                namespace numeric {
                                     namespace odesystem {


/*-----------------------------------------------------------------------------
 *       
 *    Compute Explicit (Forward) Euler solution of a given system  
 *    of ODE (RHS problem) IVP 
 *
 *    @Marco Ghiani Dec 2017, Glasgow UK
 *
 -----------------------------------------------------------------------------*/




template<typename Type = double>
class ForwardEulerSolver :
                            public Euler<Type> 
{
      
    public:  
      ForwardEulerSolver(const rhsODEProblem<Type> & that) noexcept :
                                                                        Euler<Type>{that} 
                  {}
      
      virtual ~ForwardEulerSolver() = default ;

      using OdeSolver<Type>::rhs;

      void solve(const std::string filename) override final;
      void solve() noexcept override final                 ;
//
//
  private:

      using OdeSolver<Type>::t  ;
      using OdeSolver<Type>::u  ;
      using OdeSolver<Type>::dt ; 
      using OdeSolver<Type>::t0 ;
      using OdeSolver<Type>::tf ;
      using OdeSolver<Type>::u0 ;
      
      using OdeSolver<Type>::Ns ;

};

//------------------  Implementation (to be put into .cpp file)   -----------------  //


template<typename Type>
inline void ForwardEulerSolver<Type>::solve(const std::string filename) {
      
      std::ofstream f(filename, std::ios::out );
      
      if(!f)
      {     
         throw std::runtime_error(">> Error opening File in FWD euler solver <<" );
      }
      else
      {
         std::cout << "Running ForwardEuler Solver" << std::endl;
      
         u.resize(u0().size());

         for(auto i=0; i < u0().size() ; i++ )        // set initial Value
            u[i] = u0()[i] ;
         
               
         for(t =t0() ; t <= tf() ; t+= dt() )
         {
           
           f << t << ' ' ;
           for(auto i=0 ; i < u.size() ; i++ )
              f  << u[i] <<  " " ;
           f << std::endl; 
     
           for(auto j=0 ; j< u.size() ; j++) 
               u[j] += dt() * rhs.f[j](t , u) ;
         
         } 
         std::cout << "... Done " << std::endl;  
      
         f.close();
      } 
}


template<typename Type>
inline void ForwardEulerSolver<Type>::solve() noexcept 
{
       
       
     std::cout << "Running ForwardEuler Solver" << std::endl;
      
      u.resize(u0().size());

     for(auto i=0; i < u0().size() ; i++ )      
        u[i] = u0()[i] ;                        // set initial Value 
         
               
     for(t =t0() ; t < tf() ; t+= dt() )
     {
          
       std::cout << t << ' ' ;
       for(auto i=0 ; i < u.size() ; i++ )
          std::cout  << u[i] <<  " " ;
       std::cout << std::endl; 
     
       for(auto j=0 ; j< u.size() ; j++) 
           u[j] += dt() * rhs.f[j](t , u) ;
         
     } 
     std::cout << "... Done " << std::endl;  

}
  
  }//ode
 }//numeric
}//mg
# endif
