# ifndef __CRANK_NICHOLSON_SOLVER_H__
# define __CRANK_NICHOLSON_SOLVER_H__

# include "RungeKutta.H"
# include "rhsODEproblem.H"

namespace mg {
                namespace numeric {
                                     namespace odesystem {


/*-----------------------------------------------------------------------------
 *    
 * @brief Compute Crank Nicolson Implicit scheme 
 *  using Predictor-Corrector (PC) method using explicit Euler as predictor. 
 *    
 *    (ODE) IVProblem 
 *
 *    @author Marco Ghiani 
 *    @date Dec/2017
 *    @place Glasgow UK
 * *
 *
 -----------------------------------------------------------------------------*/




template<typename Type = double>
class CrankNicholsonSolver  :
                                 public RungeKutta<Type>  
{
      
    public:  
      CrankNicholsonSolver(const rhsODEProblem<Type> & that) noexcept :
                                                                        RungeKutta<Type>{that}
                                                                 
                  {}
      
      virtual ~CrankNicholsonSolver() = default ;

      using OdeSolver<Type>::rhs;

      void solve(const std::string filename) override final ;
      void solve() noexcept override final                  ;
//
//
  private:

      using OdeSolver<Type>::t  ;
      
      using RungeKutta<Type>::u   ;
      using RungeKutta<Type>::up  ;
      using RungeKutta<Type>::uc  ;
      
      using OdeSolver<Type>::dt ; 
      using OdeSolver<Type>::t0 ;
      using OdeSolver<Type>::tf ;
      
      using OdeSolver<Type>::u0 ;
      
      using OdeSolver<Type>::Ns ;

};

//------------------  Implementation (to be put into .cpp file)   -----------------  //


template<typename Type>
inline void CrankNicholsonSolver<Type>::solve(const std::string filename) {
      
      std::ofstream f(filename, std::ios::out );
      
      if(!f)
      {     
         throw std::runtime_error(">> Error opening File in Crank Nicholson solver <<" );
      }
      else
      {
         std::cout << "Running Crank-Nicholson ( predictor - corrector ) Solver" << std::endl;
      
         u.resize (  u0().size());
         up.resize( u0().size());   
         uc.resize( u0().size());   
         for(auto i=0; i < u0().size() ; i++ )
         {
            u[i] = u0()[i] ;          // Initial Value - 
           up[i] = u0()[i] ;          // Initial Value (predictor value vector)           
           uc[i] = u0()[i] ; 
         }
                
           
         for(t=t0() ; t < tf() ; t+= dt()  )
         {
           f << t << ' ' ; 
           for(auto i=0 ; i< u.size() ; i++) 
               f << u[i] << " " ; 
           f << std::endl ;
            
           for( auto i=0 ; i < u.size() ; i++)        
           {    
              up[i] = u[i] + dt() * rhs.f[i](t , u) ;   // forward Euler PREDICTOR
           }
           for(auto i=0 ; i < u.size() ; i++)    // Crank-nicholson Corrector   
           {
               u[i] += dt()/2 * ( rhs.f[i](t , u) + rhs.f[i]( t +dt() , up ) ); 
           }             
            
         } 
         std::cout << "... Done " << std::endl;  
      
         f.close();
      } 
}

//--

template<typename Type>
inline void CrankNicholsonSolver<Type>::solve() noexcept 
{
      std::cout << "Running Crank-Nicholson ( predictor - corrector ) Solver" << std::endl;
      
         u.resize (  u0().size());
         up.resize( u0().size());   
         uc.resize( u0().size());   
         for(auto i=0; i < u0().size() ; i++ )
         {
            u[i] = u0()[i] ;          // Initial Value - 
           up[i] = u0()[i] ;          // Initial Value (predictor value vector)           
           uc[i] = u0()[i] ; 
         }
                
           
         for(t=t0() ; t < tf() ; t+= dt()  )
         {
           std::cout << t << ' ' ; 
           for(auto i=0 ; i< u.size() ; i++) 
               std::cout << u[i] << " " ; 
           std::cout << std::endl ;
            
           for( auto i=0 ; i < u.size() ; i++)        
           {      
                  up[i] = uc[i] + dt() * rhs.f[i](t , u) ;   // forward Euler PREDICTOR
                                                             // Crank-nicholson Corrector   
                  u[i] += dt()/2 * ( rhs.f[i](t , u) + rhs.f[i]( t +dt() , up ) ); 
           }             
            
         } 
         std::cout << "... Done " << std::endl;  
      
 
             
}
  
  }//ode
 }//numeric
}//mg
# endif
