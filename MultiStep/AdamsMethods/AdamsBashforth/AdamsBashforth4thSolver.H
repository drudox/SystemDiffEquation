# ifndef __ADAMS_BASHFORTH_4TH_SOLVER_H__
# define __ADAMS_BASHFORTH_4TH_SOLVER_H__

# include "../AdamsMethods.H"
# include "../../../rhsODEproblem.H"

namespace mg {
                namespace numeric {
                                     namespace odesystem {


/**-----------------------------------------------------------------------------
 *  @class AdamsBashforth4thSolver
 *
 * @brief AdamsBashforth4thSolver
 * Perform linear multi-step 4 th order (4 step) solution of a given 
 * (ODE-SDE) RHS system/problem 
 * 
 *
 * @author Marco Ghiani Dec 2017, Glasgow UK
 *
 *
 -----------------------------------------------------------------------------*/




template<typename Type = double>
class AdamsBashforth4thSolver :
                                 public AdamsMethods<Type>
{
      
    public:  

      AdamsBashforth4thSolver(const rhsODEProblem<Type> & that) noexcept :
                                                                            AdamsMethods<Type>{that} 
                  {}
      
      virtual ~AdamsBashforth4thSolver() = default ;

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
      using AdamsMethods<Type>::um2 ;
      using AdamsMethods<Type>::um3 ;
      using AdamsMethods<Type>::up1 ;

      using AdamsMethods<Type>::K1 ;
      using AdamsMethods<Type>::K2 ;
      using AdamsMethods<Type>::K3 ;
      using AdamsMethods<Type>::K4 ;
      
 

};

//------------------  Implementation (to be put into .cpp file)   -----------------  //


template<typename Type>
inline void AdamsBashforth4thSolver<Type>::solve(const std::string filename) {
      
      std::ofstream f(filename, std::ios::out );
      
      if(!f)
      {     
         throw std::runtime_error(">> Error opening File in ADAMS-BASHFORTH 4th oreder solver <<" );
      }
      else
      {
         std::cout << "Running Adams-Bashforth 4th order (4 step) Solver" << std::endl;
      
         u.resize(u0().size());
         um3.resize(u0().size());
         um2.resize(u0().size());
         um1.resize(u0().size());
         up1.resize(u0().size());
         
         K1.resize(u0().size());
         K2.resize(u0().size());
         K3.resize(u0().size());
         K4.resize(u0().size());


         t = t0();

         f << t0() << ' ' ;
         for(auto i=0 ; i < u.size() ; i++)
         {
            um3[i] = u0()[i] ;                  // set initial value 
            f << um3[i] << ' ' ;                // and write there to file
         }
         f << std::endl;
         
      
         // using for the trigger points  (start up the method)
         //
         for(auto i=0 ; i< u.size() ; i++)
         {  
            K1[i] = rhs.f[i]( t          ,  um3 );
            K2[i] = rhs.f[i]( t + dt()/2 ,  um3 + dt()/2 * K1[i] );
            K3[i] = rhs.f[i]( t + dt()/2 ,  um3 + dt()/2 * K2[i] );
            K4[i] = rhs.f[i]( t + dt()   ,  um3 + dt()   * K3[i] );
                                                 // compute first points (start - up  Adams Bashforth )
            um2[i] =  um3[i] + dt()/6 *(K1[i] + 2*K2[i] + 2*K3[i] + K4[i]) ;  // Rk 4th order 
         }
      
         t+=dt();   
            
         f << t << ' ' ;
         for(auto i=0 ; i < u.size() ; i++)
            f << um2[i] << ' ' ;                // and write there to file
         f << std::endl;
   

         for(auto i=0 ; i< u.size() ; i++)
         {  
            K1[i] = rhs.f[i]( t          ,  um2 );
            K2[i] = rhs.f[i]( t + dt()/2 ,  um2 + dt()/2*K1[i] );
            K3[i] = rhs.f[i]( t + dt()/2 ,  um2 + dt()/2*K2[i] );
            K4[i] = rhs.f[i]( t + dt()   ,  um2 + dt()  *K3[i] );
            
                                      // compute second points (start - up  Adams Bashforth )
            um1[i] =  um2[i] + dt()/6 *(K1[i] + 2*K2[i] +2*K3[i] + K4[i] ) ;  // Rk 4th order 
         }
         
         t+=dt() ;
            
         f << t << " " ;
         for(auto i=0 ; i < u.size() ; i++)
            f << um1[i] << "  " ;
         f << std::endl;   


         for(auto i=0; i < u.size() ; i++)
         {
            K1[i] = rhs.f[i]( t            , um1 );
            K2[i] = rhs.f[i]( t +   dt()/2 , um1+ dt()/2 *K1[i] );
            K3[i] = rhs.f[i]( t +   dt()/2 , um1+ dt()/2 *K2[i] );
            K4[i] = rhs.f[i]( t +   dt()   , um1+ dt()   *K3[i] );
                                                             // compute third points
            u[i]   = um1[i] + dt()/6 *(K1[i] + 2*K2[i] + 2*K3[i] + K4[i]);  // RK 4th order 
         }

         for(t +=dt() ; t <= tf() ; t += dt() )
         {
            f << t << "  " ;
            for(auto i=0 ; i < u.size() ; i++)
                  f << u[i] << " " ;
            f << std::endl; 
            
            for(auto i=0 ; i< u.size(); i++)
            { 
               up1[i] = u[i] + dt()/24.0 *( 55.0 * rhs.f[i](t , u) 
                                          - 59.0 * rhs.f[i](t-dt() , um1) 
                                          + 37.0 * rhs.f[i](t-(2*dt()), um2 ) 
                                          -  9.0 * rhs.f[i](t-(3*dt()), um3 ) );
 
                  
            }
            for(auto i=0 ; i < u.size() ; i++)
            {
               um3[i] = um2[i]; 
               um2[i] = um1[i];     
               um1[i] =   u[i];
                 u[i] = up1[i];
            }
          
         }  
         std::cout << "... Done " << std::endl;  
         f.close();
      } 
}


template<typename Type>
inline void AdamsBashforth4thSolver<Type>::solve() noexcept 
{
    std::cout << "Running Adams-Bashforth 4th order (4 step) Solver" << std::endl;
      
         u.resize(u0().size());
         um3.resize(u0().size());
         um2.resize(u0().size());
         um1.resize(u0().size());
         up1.resize(u0().size());
         
         K1.resize(u0().size());
         K2.resize(u0().size());
         K3.resize(u0().size());
         K4.resize(u0().size());


         t = t0();

         std::cout << t0() << ' ' ;
         for(auto i=0 ; i < u.size() ; i++)
         {
            um3[i] = u0()[i] ;                  // set initial value 
            std::cout << um2[i] << ' ' ;                // and write there to file
         }
         std::cout << std::endl;
         
      
         // using for the trigger points  (start up the method)
         //
         for(auto i=0 ; i< u.size() ; i++)
         {  
            K1[i] = rhs.f[i]( t          ,  um3 );
            K2[i] = rhs.f[i]( t + dt()/2 ,  um3 + dt()/2 * K1[i] );
            K3[i] = rhs.f[i]( t + dt()/2 ,  um3 + dt()/2 * K2[i] );
            K4[i] = rhs.f[i]( t + dt()   ,  um3 + dt()   * K3[i] );
                                                 // compute first points (start - up  Adams Bashforth )
            um2[i] =  um3[i] + dt()/6 *(K1[i] + 2*K2[i] + 2*K3[i] + K4[i]) ;  // Rk 4th order 
         }
      
         t+=dt();   
            
         std::cout << t << ' ' ;
         for(auto i=0 ; i < u.size() ; i++)
            std::cout << um2[i] << ' ' ;                // and write there to file
         std::cout << std::endl;
   

         for(auto i=0 ; i< u.size() ; i++)
         {  
            K1[i] = rhs.f[i]( t          ,  um2 );
            K2[i] = rhs.f[i]( t + dt()/2 ,  um2 + dt()/2*K1[i] );
            K3[i] = rhs.f[i]( t + dt()/2 ,  um2 + dt()/2*K2[i] );
            K4[i] = rhs.f[i]( t + dt()   ,  um2 + dt()  *K3[i] );
            
                                      // compute second points (start - up  Adams Bashforth )
            um1[i] =  um2[i] + dt()/6 *(K1[i] + 2*K2[i] +2*K3[i] + K4[i] ) ;  // Rk 4th order 
         }
         
         t+=dt() ;
            
         std::cout << t << " " ;
         for(auto i=0 ; i < u.size() ; i++)
            std::cout << um1[i] << "  " ;
         std::cout << std::endl;   


         for(auto i=0; i < u.size() ; i++)
         {
            K1[i] = rhs.f[i]( t            , um1 );
            K2[i] = rhs.f[i]( t +   dt()/2 , um1+ dt()/2 *K1[i] );
            K3[i] = rhs.f[i]( t +   dt()/2 , um1+ dt()/2 *K2[i] );
            K4[i] = rhs.f[i]( t +   dt()   , um1+ dt()   *K3[i] );
                                                             // compute third points
            u[i]   = um1[i] + dt()/6 *(K1[i] + 2*K2[i] + 2*K3[i] + K4[i]);  // RK 4th order 
         }

         for(t +=dt() ; t <= tf() ; t += dt() )
         {
            std::cout << t << "  " ;
            for(auto i=0 ; i < u.size() ; i++)
                  std::cout << u[i] << " " ;
            std::cout << std::endl; 
            
            for(auto i=0 ; i< u.size(); i++)
            { 
               up1[i] = u[i] + dt()/24.0 *( 55.0 * rhs.f[i](t , u) 
                                          - 59.0 * rhs.f[i](t-dt() , um1) 
                                          + 37.0 * rhs.f[i](t-(2*dt()), um2 ) 
                                          -  9.0 * rhs.f[i](t-(3*dt()), um3 ) );
 
                  
            }
            for(auto i=0 ; i < u.size() ; i++)
            {
               um3[i] = um2[i]; 
               um2[i] = um1[i];     
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
