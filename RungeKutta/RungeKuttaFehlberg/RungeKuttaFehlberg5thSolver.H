# ifndef __RUNGEKUTTA_FEHLBERG_5ORD_SOLVER_H__
# define __RUNGEKUTTA_FEHLBERG_5ORD_SOLVER_H__

# include "rhsODEproblem.H"
# include "RungeKutta.H" 

namespace mg { 
                namespace numeric {
                                    namespace odesystem {


/*-------------------------------------------------------------------------------
 *    
 *    @class RKFehlberg54thSolver :
 *
 *    Perform Numerical solution of a System, of first order ODE du/dt = f(u,t)
 *    this is an explicit Runge Kutta Method of 5 order due to Dormand & Prince
 *    
 *    (whit step-size control) 
 *    
 *
 *    @author Marco Ghiani Dec 2017, Glasgow UK
 *
 ------------------------------------------------------------------------------*/



template<typename Type= double>
class RKFehlberg54thSolver              // Dormand Prince Solver (explicit Runge Kutta )
                         :   public  RungeKutta<Type>    
{
      
    public:  
      RKFehlberg54thSolver(const rhsODEProblem<Type> & that) noexcept :
                                                                        RungeKutta<Type>{that},
                                                                        dt0{that.dt()}
                  {}
      
      virtual ~RKFehlberg54thSolver() = default;
      

      using OdeSolver<Type>::rhs;

      void solve(const std::string filename) override final ;
      void solve() noexcept override final                     ;
     
    private:
      
      using OdeSolver<Type>::t  ;
      using OdeSolver<Type>::u  ;
      using OdeSolver<Type>::dt ; 
      using OdeSolver<Type>::t0 ;
      using OdeSolver<Type>::tf ;
      using OdeSolver<Type>::u0 ;
      
      using OdeSolver<Type>::Ns ;
      
      using RungeKutta<Type>::u_4 ;
      using RungeKutta<Type>::u_5 ;

      using OdeSolver<Type>::toll ;
      using OdeSolver<Type>::setTimeStep ;

      using RungeKutta<Type>::k1;
      using RungeKutta<Type>::k2;
      using RungeKutta<Type>::k3;
      using RungeKutta<Type>::k4;
      using RungeKutta<Type>::k5;
      using RungeKutta<Type>::k6;
 
      using RungeKutta<Type>::s;
      
      const Type dt0;
      
      using RungeKutta<Type>::stepToll;
      
      using RungeKutta<Type>::truncError;

};


template<typename Type>
inline void RKFehlberg54thSolver<Type>::solve(const std::string filename)  {
      
      std::ofstream f(filename, std::ios::out );
      
      if(!f)
      {     
         std::string mess = "Error opening file " + filename + " in R-K Fehlberg 5th order Solver " ;  
         throw std::runtime_error(mess.c_str());   
      }
      else
      {
         std::cout << "Running Runge-Kutta Felhberg 5th order Solver" << std::endl;
         

         t = t0();
         u.resize(u0().size());
         u_4.resize(u0().size());
         u_5.resize(u0().size());
         
         truncError.resize(u0().size());
         s.resize(u0().size());
         
         k1.resize(u0().size());
         k2.resize(u0().size());
         k3.resize(u0().size());
         k4.resize(u0().size());
         k5.resize(u0().size());
         k6.resize(u0().size());

        // std::vector<Type> Dt(u0().size(), dt());
         
         for(auto i=0 ; i < u0().size() ; i++)
         {
            u[i]   = u0()[i];                           //set Init Value  
            u_5[i] = u0()[i];
         }
         for(t = t0() ;  t <= tf() ; t += dt() )
         {
            f << t << ' ' ;
            for(auto i=0 ; i< u.size() ; i++)                  // write Initial Value 
               f << u[i] << " "; 
            f << std::endl;
            
           for(auto j=0 ; j< u.size() ; j++)
           { 
            
             k1[j] = dt()*rhs.f[j](t                  , u_5                                             );   
             k2[j] = dt()*rhs.f[j](t+dt()/4           , u_5 + k1[j]/4.0                                 ); 
             k3[j] = dt()*rhs.f[j](t+ 3.0/8.0*dt()    , u_5 + 3.0/32.0*k1[j] + 9.0/32.0* k2[j]          );
             k4[j] = dt()*rhs.f[j](t+ 12.0/13.0* dt() , u_5 + 1932.0/2197.0* k1[j] -  7200.0/2197.0*k2[j] 
                                                            + 7296.0/2197.0* k3[j]                      ); 
             k5[j] = dt()*rhs.f[j](t+dt()             , u_5 + 439.0/216.0  * k1[j] - 8.0*k2[j]  
                                                            + 3680.0/513.0 * k3[j] - 845.0/4104.0*k4[j] ); 
             k6[j] = dt()*rhs.f[j](t+dt()/2.0         , u_5 - 8.0/27.0 * k1[j] + 2.0*k2[j] 
                                                            - 3544.0/2565.0* k3[j] + 1859.0/4104.0*k4[j]
                                                                                   - 11.0/40.0    *k5[j]); 
                  

            // 5-th order Runge Kutta (Fehlberg ) method 
             u_5[j] = u[j] + (16.0/135.0*k1[j] + 6656.0/12825.0 *k3[j] + 28561.0/56430.0* k4[j] 
                                                 - 9.0/50.0* k5[j] + 2.0/55.0* k6[j]) ;  
            
            // 4-th order Runge Kutta (Fehlberg) method
             u_4[j] = u[j] + (25.0/216.0*k1[j] + 1408.0/2565.0*k3[j] + 2197.0/4104.0*k4[j]-1.0/5.0*k5[j]);   

             truncError[j] = fabs(u_5[j]-u_4[j]);
            
             s[j] = 0.8409*(pow(stepToll*dt()/truncError[j],(1./4.)));  // adaptative Time-step
                 
             u[j] = u_5[j] ;
                  
            }
            /*
            for(auto j=0; j< u.size() ; j++)
            { 
              if( (s.min()*dt() ) < dt0 ){
                setTimeStep( s.min() * dt()/2 );  
                std::cout << "t1  " <<s.min()*dt() << std::endl; // adaptative time step
              }
              else if( s.min()*dt() > dt0 ){
                 setTimeStep(dt0/2);
                 std::cout << dt0 << std::endl; // adaptative time step
              }
            }
            */
         } 
         std::cout << "... Done " << std::endl;  
      
      f.close();
      }
}


template<typename Type>
inline void RKFehlberg54thSolver<Type>::solve() noexcept 
{
        std::cout << "Running Runge-Kutta Felhberg 5th order Solver" << std::endl;
         

         t = t0();
         u.resize(u0().size());
         k1.resize(u0().size());
         k2.resize(u0().size());
         k3.resize(u0().size());
         k4.resize(u0().size());
         k5.resize(u0().size());
         k6.resize(u0().size());

        // std::vector<Type> Dt(u0().size(), dt());
         
         for(auto i=0 ; i < u0().size() ; i++)
            u[i] = u0()[i];                           //set Init Value  
         
         for(t = t0() ;  t <= tf() ; t += dt() )
         {
            std::cout << t << ' ' ;
            for(auto i=0 ; i< u.size() ; i++)                  // write Initial Value 
               std::cout << u[i] << " "; 
            std::cout  << std::endl;

           for(auto j=0 ; j< u.size() ; j++)
           { 
            
             k1[j] = dt()*rhs.f[j](t                  , u                                               );   
             k2[j] = dt()*rhs.f[j](t+dt()/4           , u + k1[j]/4.0                                   ); 
             k3[j] = dt()*rhs.f[j](t+ 3.0/8.0*dt()    , u + 3.0/32.0*k1[j] + 9.0/32.0* k2[j]            );
             k4[j] = dt()*rhs.f[j](t+ 12.0/13.0* dt() , u + 1932.0/2197.0* k1[j] -  7200.0/2197.0*k2[j] 
                                                          + 7296.0/2197.0* k3[j]                        ); 
             k5[j] = dt()*rhs.f[j](t+dt()             , u + 439.0/216.0  * k1[j] - 8.0*k2[j]  
                                                          + 3680.0/513.0 * k3[j] - 845.0/4104.0*k4[j]   ); 
             k6[j] = dt()*rhs.f[j](t+dt()/2.0         , u - 8.0/27.0 * k1[j] + 2.0*k2[j] 
                                                          - 3544.0/2565.0* k3[j] + 1859.0/4104.0*k4[j]
                                                                                 - 11.0/40.0    *k5[j]  ); 
                  

            // 5-th order Runge Kutta (Fehlberg ) method 
              u[j] += (16.0/135.0*k1[j] + 6656.0/12825.0 *k3[j] + 28561.0/56430.0* k4[j] - 9.0/50.0* k5[j] 
                                                                                     + 2.0/55.0* k6[j]) ;  
            }
           /*
            for(auto j=0; j< u.size() ; j++)
            { 
              if( (s.min()*dt() ) < dt0 ){
                setTimeStep( s.min() * dt()/2 );  
                std::cout << "t1  " <<s.min()*dt() << std::endl; // adaptative time step
              }
              else if( s.min()*dt() > dt0 ){
                 setTimeStep(dt0/2);
                 std::cout << dt0 << std::endl; // adaptative time step
              }
            }
            */
 
         } 
         std::cout << "... Done " << std::endl;  

}

  }//ode
 }//numeric
}//mg 
# endif
