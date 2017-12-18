# ifndef __ADAMS_BASHFORTH_2ND_SOLVER_H__
# define __ADAMS_BASHFORTH_2ND_SOLVER_H__

# include "../AdamsMethods.H"
# include "../../../rhsODEproblem.H"

namespace mg {
                namespace numeric {
                                     namespace odesystem {


/**-----------------------------------------------------------------------------
 *
 *\ @brief \class AdamsBashforth2ndSolver
 *\ Perform multi-step 2 order accuracy solution of a given (ODE) RHS problem 
 * 
 *
 *\ @Marco Ghiani Dec 2017, Glasgow UK
 *
 *
 -----------------------------------------------------------------------------*/




template<typename Type = double>
class AdamsBashforth2ndSolver :
                                 public AdamsMethods<Type> 
{
      
    public:  
      AdamsBashforth2ndSolver(const rhsODEProblem<Type> & that) noexcept :
                                                                           AdamsMethods<Type>{that} 
                  {}
      
      virtual ~AdamsBashforth2ndSolver() = default ;

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
      
      using AdamsMethods<Type>::um1 ;
      using AdamsMethods<Type>::up1 ;
      
      using AdamsMethods<Type>::K1 ;
      using AdamsMethods<Type>::K2 ;


};

//------------------  Implementation (to be put into .cpp file)   -----------------  //


template<typename Type>
inline void AdamsBashforth2ndSolver<Type>::solve(const std::string filename) {
      
      std::ofstream f(filename, std::ios::out );
      
      if(!f)
      {     
         throw std::runtime_error(">> Error opening File in ADAMS-BASHFORTH 2step solver <<" );
      }
      else
      {
         std::cout << "Running Adams-Bashforth 2nd order (2 step) Solver" << std::endl;
      
         u.resize(u0().size());
         um1.resize(u0().size());
         up1.resize(u0().size());
         
         K1.resize(u.size());
         K2.resize(u.size());
           
         t = t0();

         f << t << ' ' ;
         for(auto i=0 ; i < u.size() ; i++)
         {
            um1[i] = u0()[i] ;                  // set initial value 
            f << um1[i] << ' ' ;                // and write there to file
         }
         f << std::endl;
         

         for(auto i=0 ; i< u.size() ; i++)
         {  
            K1[i] = rhs.f[i]( t        ,  um1 );
            K2[i] = rhs.f[i]( t + dt() ,  um1 + dt()*K1[i] );
                                      // compute first point (start - up  Adams Bashforth )
            u[i] =  um1[i] + dt()/2 *(K1[i] + K2[i]) ;  // Rk 2nd order 
         }

        /*
         for(auto i=0 ; i< u.size() ; i++)
         {  
            k1[i] = rhs.f[i]( t          ,  um1 );
            k2[i] = rhs.f[i]( t + dt()/2 ,  um1 + k1*dt()/2 );
            
            // initiation first point 
            u[i] =  um1[i] + dt() * k2[i];  // Rk 2nd order PREDICTOR
         }
      */
         for(t=t0()+dt() ; t <= tf() ; t += dt() )
         {
            f << t << "  " ;
            for(auto i=0 ; i < u.size() ; i++)
                  f << u[i] << " " ;
            f << std::endl; 
            
            for(auto i=0 ; i< u.size(); i++)
            { 
                  up1[i] = u[i] + dt()/2 *(3 * rhs.f[i](t , u) - rhs.f[i](t-dt() , um1))  ;
            }
            for(auto i=0 ; i< u.size(); i++)
            {
               um1[i] =   u[i];
                 u[i] = up1[i];
            } 
         }  
         std::cout << "... Done " << std::endl;  
         f.close();
      } 
}


template<typename Type>
inline void AdamsBashforth2ndSolver<Type>::solve() noexcept 
{
        std::cout << "Running Adams-Bashforth (2 step) Solver" << std::endl;
      
         u.resize(u0().size());
         um1.resize(u0().size());
         up1.resize(u0().size());
         

         std::cout << t0() << ' ' ;
         for(auto i=0 ; i < u.size() ; i++)
         {
            um1[i] = u0()[i] ;                  // set initial value 
            std::cout << um1[i] << ' ' ;                // and write there to file
         }
         std::cout << std::endl;
         


         std::valarray<Type> k1(u.size());
         std::valarray<Type> k2(u.size());
         
         for(auto i=0 ; i< u.size() ; i++)
         {  
            k1[i] = rhs.f[i]( t          ,  um1 );
            k2[i] = rhs.f[i]( t + dt()/2 ,  um1 + k1*dt()/2 );
            
            // initiation first point 
            u[i] =  um1[i] + dt() * k2[i];  // Rk 2nd order PREDICTOR
         }

         for(t=t0()+dt() ; t < tf() ; t += dt() )
         {
            std::cout << t << "  " ;
            for(auto i=0 ; i < u.size() ; i++)
                  std::cout << u[i] << " " ;
            std::cout << std::endl; 
            
            for(auto i=0 ; i< u.size(); i++)
            { 
               up1[i] = u[i] + dt()/2 *(3 * rhs.f[i](t , u) - rhs.f[i](t-dt() , um1))  ;
               
               um1[i] =   u[i];
                 u[i] = up1[i];
            } 
         }  
         std::cout << "... Done " << std::endl;  
 
}
  
  }//ode
 }//numeric
}//mg
# endif
