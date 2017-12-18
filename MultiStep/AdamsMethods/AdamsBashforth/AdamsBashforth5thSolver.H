# ifndef __ADAMS_BASHFORTH_5TH_SOLVER_H__
# define __ADAMS_BASHFORTH_5TH_SOLVER_H__

# include "../AdamsMethods.H"
# include "../../../rhsODEproblem.H"

namespace mg {
                namespace numeric {
                                     namespace odesystem {


/**-----------------------------------------------------------------------------
 *  @class AdamsBashforth5thSolver
 *
 * @brief AdamsBashforth5thSolver
 * Perform linear multi-step 4 th order (4 step) solution of a given 
 * (ODE-SDE) RHS system/problem 
 * 
 *
 * @author Marco Ghiani Dec 2017, Glasgow UK
 *
 *
 -----------------------------------------------------------------------------*/




template<typename Type = double>
class AdamsBashforth5thSolver :
                                 public AdamsMethods<Type>
{
      
    public:  

      AdamsBashforth5thSolver(const rhsODEProblem<Type> & that) noexcept :
                                                                            AdamsMethods<Type>{that} 
                  {}
      
      virtual ~AdamsBashforth5thSolver() = default ;

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
      using AdamsMethods<Type>::um4 ;
      using AdamsMethods<Type>::up1 ;

      using AdamsMethods<Type>::K1 ;
      using AdamsMethods<Type>::K2 ;
      using AdamsMethods<Type>::K3 ;
      using AdamsMethods<Type>::K4 ;
      using AdamsMethods<Type>::K5 ;
      
 

};

//------------------  Implementation (to be put into .cpp file)   -----------------  //


template<typename Type>
inline void AdamsBashforth5thSolver<Type>::solve(const std::string filename) {
      
      std::ofstream f(filename, std::ios::out );
      
      if(!f)
      {     
         throw std::runtime_error(">> Error opening File in ADAMS-BASHFORTH 5th oreder solver <<" );
      }
      else
      {
         std::cout << "Running Adams-Bashforth 5th order (5 step) Solver" << std::endl;
      
         u.resize(u0().size());
         um4.resize(u0().size());
         um3.resize(u0().size());
         um2.resize(u0().size());
         um1.resize(u0().size());
         up1.resize(u0().size());
         
         K1.resize(u0().size());
         K2.resize(u0().size());
         K3.resize(u0().size());
         K4.resize(u0().size());
         K5.resize(u0().size());


         t = t0();

         f << t0() << ' ' ;
         for(auto i=0 ; i < u.size() ; i++)
         {
            um4[i] = u0()[i] ;                  // set initial value 
            f << um4[i] << ' ' ;                // and write there to file
         }
         f << std::endl;
         
      
         for(auto i=0 ; i< u.size() ; i++)
         {  
            K1[i] = dt()*rhs.f[i]( t            ,  um4                                      );
            K2[i] = dt()*rhs.f[i]( t + dt()/3.0 ,  um4 + K1[i]/3.0                          );
            K3[i] = dt()*rhs.f[i]( t + dt()/3.0 ,  um4 + 1.0/6.0*(K1[i] + K2[i])            );
            K4[i] = dt()*rhs.f[i]( t + dt()/2.0 ,  um4 + 1.0/8.0*(K1[i] + 3.0*K3[i])        );
            K5[i] = dt()*rhs.f[i]( t + dt()     ,  um4 + 1.0/2.0*(K1[i] - 3.0*K3[i]+ 4.0*K4));

                                                      // compute first points (start - up  Adams Bashforth )
            um3[i] =  um4[i] + 1.0/6.0 *(K1[i] + 4*K4[i] + 1.0*K5[i] ) ; // Runge-Kutta-Merson 5th order 
         }
   
         t+=dt();   
            
         f << t << ' ' ;
         for(auto i=0 ; i < u.size() ; i++)
            f << um3[i] << ' ' ;                // write them to file
         f << std::endl;
 
         



         // using for the trigger points  (start up the method)
         //
         for(auto i=0 ; i< u.size() ; i++)
         {  
            K1[i] = dt()*rhs.f[i]( t            ,  um3                                      );
            K2[i] = dt()*rhs.f[i]( t + dt()/3.0 ,  um3 + K1[i]/3.0                          );
            K3[i] = dt()*rhs.f[i]( t + dt()/3.0 ,  um3 + 1.0/6.0*(K1[i] + K2[i])            );
            K4[i] = dt()*rhs.f[i]( t + dt()/2.0 ,  um3 + 1.0/8.0*(K1[i] + 3.0*K3[i])        );
            K5[i] = dt()*rhs.f[i]( t + dt()     ,  um3 + 1.0/2.0*(K1[i] - 3.0*K3[i]+ 4.0*K4));

                                            /// compute second points (start - up  Adams Bashforth )
            um2[i] =  um3[i] + 1.0/6.0 *(K1[i] + 4*K4[i] + 1.0*K5[i] ) ; // Runge-Kutta-Merson 5th order 
         }
      
         t+=dt();   
            
         f << t << ' ' ;
         for(auto i=0 ; i < u.size() ; i++)
            f << um2[i] << ' ' ;                // and write there to file
         f << std::endl;
   

         for(auto i=0 ; i< u.size() ; i++)
         {  
            K1[i] = dt()*rhs.f[i]( t            ,  um2                                      );
            K2[i] = dt()*rhs.f[i]( t + dt()/3.0 ,  um2 + K1[i]/3.0                          );
            K3[i] = dt()*rhs.f[i]( t + dt()/3.0 ,  um2 + 1.0/6.0*(K1[i] + K2[i])            );
            K4[i] = dt()*rhs.f[i]( t + dt()/2.0 ,  um2 + 1.0/8.0*(K1[i] + 3.0*K3[i])        );
            K5[i] = dt()*rhs.f[i]( t + dt()     ,  um2 + 1.0/2.0*(K1[i] - 3.0*K3[i]+ 4.0*K4));

                                      /// compute third points (start - up  Adams Bashforth )
            um1[i] =  um2[i] + 1.0/6.0 *(K1[i] + 4*K4[i] + 1.0*K5[i] ) ; // Runge-Kutta-Merson 5th order 
         }
         
         t+=dt() ;
            
         f << t << " " ;
         for(auto i=0 ; i < u.size() ; i++)
            f << um1[i] << "  " ;
         f << std::endl;   


         for(auto i=0; i < u.size() ; i++)
         {
            K1[i] = dt()*rhs.f[i]( t            ,  um1                                      );
            K2[i] = dt()*rhs.f[i]( t + dt()/3.0 ,  um1 + K1[i]/3.0                          );
            K3[i] = dt()*rhs.f[i]( t + dt()/3.0 ,  um1 + 1.0/6.0*(K1[i] + K2[i])            );
            K4[i] = dt()*rhs.f[i]( t + dt()/2.0 ,  um1 + 1.0/8.0*(K1[i] + 3.0*K3[i])        );
            K5[i] = dt()*rhs.f[i]( t + dt()     ,  um1 + 1.0/2.0*(K1[i] - 3.0*K3[i]+ 4.0*K4));

                                      /// compute fourth points (start - up  Adams Bashforth )
            u[i] =  um1[i] + 1.0/6.0 *(K1[i] + 4*K4[i] + 1.0*K5[i] ) ; // Runge-Kutta-Merson 5th order 
         }


      ///@ Main LOOP
         for(t +=dt() ; t <= tf() ; t += dt() )
         {
            f << t << "  " ;
            for(auto i=0 ; i < u.size() ; i++)
                  f << u[i] << " " ;
            f << std::endl; 
            
            for(auto i=0 ; i< u.size(); i++)
            { 
               up1[i] = u[i] + dt()       *( 1901.0/720.0 * rhs.f[i](t , u) 
                                            -1387.0/360.0 * rhs.f[i](t-dt() , um1) 
                                            +109.0/30.0 * rhs.f[i](t-(2*dt()), um2 )  
                                            -637.0/360.0 * rhs.f[i](t-(3*dt()), um3 ) 
                                             +251.0/720.0  * rhs.f[i](t-(4*dt()), um4 )  );
            }
            for(auto i=0 ; i < u.size() ; i++)
            {
               um4[i] = um3[i];  
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
inline void AdamsBashforth5thSolver<Type>::solve() noexcept 
{
    std::cout << "Running Adams-Bashforth 5th order (5 step) Solver" << std::endl;
      
         u.resize(u0().size());
         um4.resize(u0().size());
         um3.resize(u0().size());
         um2.resize(u0().size());
         um1.resize(u0().size());
         up1.resize(u0().size());
         
         K1.resize(u0().size());
         K2.resize(u0().size());
         K3.resize(u0().size());
         K4.resize(u0().size());
         K5.resize(u0().size());


         t = t0();

         std::cout << t0() << ' ' ;
         for(auto i=0 ; i < u.size() ; i++)
         {
            um4[i] = u0()[i] ;                  // set initial value 
            std::cout << um4[i] << ' ' ;                // and write there to file
         }
         std::cout << std::endl;
         
      
         for(auto i=0 ; i< u.size() ; i++)
         {  
            K1[i] = dt()*rhs.f[i]( t            ,  um4                                      );
            K2[i] = dt()*rhs.f[i]( t + dt()/3.0 ,  um4 + K1[i]/3.0                          );
            K3[i] = dt()*rhs.f[i]( t + dt()/3.0 ,  um4 + 1.0/6.0*(K1[i] + K2[i])            );
            K4[i] = dt()*rhs.f[i]( t + dt()/2.0 ,  um4 + 1.0/8.0*(K1[i] + 3.0*K3[i])        );
            K5[i] = dt()*rhs.f[i]( t + dt()     ,  um4 + 1.0/2.0*(K1[i] - 3.0*K3[i]+ 4.0*K4));

                                                      // compute first points (start - up  Adams Bashforth )
            um3[i] =  um4[i] + 1.0/6.0 *(K1[i] + 4*K4[i] + 1.0*K5[i] ) ; // Runge-Kutta-Merson 5th order 
         }
   
         t+=dt();   
            
         std::cout << t << ' ' ;
         for(auto i=0 ; i < u.size() ; i++)
            std::cout << um3[i] << ' ' ;                // write them to file
         std::cout << std::endl;
 
         



         // using for the trigger points  (start up the method)
         //
         for(auto i=0 ; i< u.size() ; i++)
         {  
            K1[i] = dt()*rhs.f[i]( t            ,  um3                                      );
            K2[i] = dt()*rhs.f[i]( t + dt()/3.0 ,  um3 + K1[i]/3.0                          );
            K3[i] = dt()*rhs.f[i]( t + dt()/3.0 ,  um3 + 1.0/6.0*(K1[i] + K2[i])            );
            K4[i] = dt()*rhs.f[i]( t + dt()/2.0 ,  um3 + 1.0/8.0*(K1[i] + 3.0*K3[i])        );
            K5[i] = dt()*rhs.f[i]( t + dt()     ,  um3 + 1.0/2.0*(K1[i] - 3.0*K3[i]+ 4.0*K4));

                                            /// compute second points (start - up  Adams Bashforth )
            um2[i] =  um3[i] + 1.0/6.0 *(K1[i] + 4*K4[i] + 1.0*K5[i] ) ; // Runge-Kutta-Merson 5th order 
         }
      
         t+=dt();   
            
         std::cout << t << ' ' ;
         for(auto i=0 ; i < u.size() ; i++)
            std::cout << um2[i] << ' ' ;                // and write there to file
         std::cout << std::endl;
   

         for(auto i=0 ; i< u.size() ; i++)
         {  
            K1[i] = dt()*rhs.f[i]( t            ,  um2                                      );
            K2[i] = dt()*rhs.f[i]( t + dt()/3.0 ,  um2 + K1[i]/3.0                          );
            K3[i] = dt()*rhs.f[i]( t + dt()/3.0 ,  um2 + 1.0/6.0*(K1[i] + K2[i])            );
            K4[i] = dt()*rhs.f[i]( t + dt()/2.0 ,  um2 + 1.0/8.0*(K1[i] + 3.0*K3[i])        );
            K5[i] = dt()*rhs.f[i]( t + dt()     ,  um2 + 1.0/2.0*(K1[i] - 3.0*K3[i]+ 4.0*K4));

                                      /// compute third points (start - up  Adams Bashforth )
            um1[i] =  um2[i] + 1.0/6.0 *(K1[i] + 4*K4[i] + 1.0*K5[i] ) ; // Runge-Kutta-Merson 5th order 
         }
         
         t+=dt() ;
            
         std::cout << t << " " ;
         for(auto i=0 ; i < u.size() ; i++)
            std::cout << um1[i] << "  " ;
         std::cout << std::endl;   


         for(auto i=0; i < u.size() ; i++)
         {
            K1[i] = dt()*rhs.f[i]( t            ,  um1                                      );
            K2[i] = dt()*rhs.f[i]( t + dt()/3.0 ,  um1 + K1[i]/3.0                          );
            K3[i] = dt()*rhs.f[i]( t + dt()/3.0 ,  um1 + 1.0/6.0*(K1[i] + K2[i])            );
            K4[i] = dt()*rhs.f[i]( t + dt()/2.0 ,  um1 + 1.0/8.0*(K1[i] + 3.0*K3[i])        );
            K5[i] = dt()*rhs.f[i]( t + dt()     ,  um1 + 1.0/2.0*(K1[i] - 3.0*K3[i]+ 4.0*K4));

                                      /// compute fourth points (start - up  Adams Bashforth )
            u[i] =  um1[i] + 1.0/6.0 *(K1[i] + 4*K4[i] + 1.0*K5[i] ) ; // Runge-Kutta-Merson 5th order 
         }


      ///@ Main LOOP
         for(t +=dt() ; t <= tf() ; t += dt() )
         {
            std::cout << t << "  " ;
            for(auto i=0 ; i < u.size() ; i++)
                  std::cout << u[i] << " " ;
            std::cout << std::endl; 
            
            for(auto i=0 ; i< u.size(); i++)
            { 
               up1[i] = u[i] + dt()       *( 1901.0/720.0 * rhs.f[i](t         , u   ) 
                                            -1387.0/360.0 * rhs.f[i](t-dt()    , um1 ) 
                                            +109.0/30.0   * rhs.f[i](t-(2*dt()), um2 )  
                                            -637.0/360.0  * rhs.f[i](t-(3*dt()), um3 ) 
                                            +251.0/720.0  * rhs.f[i](t-(4*dt()), um4 )  );

            }
            for(auto i=0 ; i < u.size() ; i++)
            {
               um4[i] = um3[i];   
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
