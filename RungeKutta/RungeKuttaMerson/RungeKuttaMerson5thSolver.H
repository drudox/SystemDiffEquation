# ifndef __RUNGEKUTTA_MERSON_5ORD_SOLVER_H__
# define __RUNGEKUTTA_MERSON_5ORD_SOLVER_H__

//# include "OdeSolver.H"
# include "rhsODEproblem.H"
# include "RungeKutta.H" 

namespace mg { 
                namespace numeric {
                                    namespace odesystem {


/*-------------------------------------------------------------------------------
 *    
 *    @class RungeKuttaMerson5thSolver :
 *    
 *    Perform Numerical solution of a System, of first order ODE du/dt = f(u,t)
 *    this is an explicit Runge Kutta Method of 5 order due to Merson
 *    (whitout step-size control) 
 *    
 *
 *    @author Marco Ghiani Dec 2017, Glasgow UK
 *
 ------------------------------------------------------------------------------*/



template<typename Type= double>
class RungeKuttaMerson5thSolver              // Dormand Prince Solver (explicit Runge Kutta )
                         :   public  RungeKutta<Type>    
{
      
    public:  
      RungeKuttaMerson5thSolver(const rhsODEProblem<Type> & that) noexcept :
                                                                        RungeKutta<Type>{that} 
                  {}
      
      virtual ~RungeKuttaMerson5thSolver() = default;
      

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
      
      using OdeSolver<Type>::toll ;

      using RungeKutta<Type>::k1;
      using RungeKutta<Type>::k2;
      using RungeKutta<Type>::k3;
      using RungeKutta<Type>::k4;
      using RungeKutta<Type>::k5;
 

};


template<typename Type>
inline void RungeKuttaMerson5thSolver<Type>::solve(const std::string filename)  {
      
      std::ofstream f(filename, std::ios::out );
      
      if(!f)
      {     
         std::string mess = "Error opening file " + filename + " in Runge-Kutta-Merson 5th order Solver " ;  
         throw std::runtime_error(mess.c_str());   
      }
      else
      {
         std::cout << "Running Runge-Kutta-Merson 5th order Solver" << std::endl;
         

         t = t0();
         u.resize(u0().size());
         k1.resize(u0().size());
         k2.resize(u0().size());
         k3.resize(u0().size());
         k4.resize(u0().size());
         k5.resize(u0().size());

        // std::vector<Type> Dt(u0().size(), dt());
         
         for(auto i=0 ; i < u0().size() ; i++)
            u[i] = u0()[i];                           //set Init Value  
         
         for(t = t0() ;  t <= tf() ; t += dt() )
         {
            f << t << ' ' ;
            for(auto i=0 ; i< u.size() ; i++)                  // write Initial Value 
               f << u[i] << " "; 
            f << std::endl;

           for(auto j=0 ; j< u.size() ; j++)
           { 
            
             // dt() =   <- TO DO   

             k1[j] = dt()*rhs.f[j](t                  , u                                          );   
             k2[j] = dt()*rhs.f[j](t+ 1.0/3.0*dt()    , u + k1[j]/3.0                              ); 
             k3[j] = dt()*rhs.f[j](t+ 1.0/3.0*dt()    , u + 1.0/6.0*(k1[j] + k2[j])                );
             k4[j] = dt()*rhs.f[j](t+ 1.0/2.0*dt()    , u + 1.0/8.0*(k1[j] + 3.0*k3[j] )           ); 
             k5[j] = dt()*rhs.f[j](t+dt()             , u + 1.0/2.0*(k1[j] - 3.0*k3[j] + 4.0*k4[j])); 
                  

            // 5-th order Runge Kutta (Fehlberg ) method 
              u[j] += 1.0/6.0*( k1[j] + 4.0*k4[j] + 1.0*k5[j] ) ;  
            }
            
         } 
         std::cout << "... Done " << std::endl;  
      
      f.close();
      }
}


template<typename Type>
inline void RungeKuttaMerson5thSolver<Type>::solve() noexcept 
{
      std::cout << "Running Runge-Kutta-Merson 5th order Solver" << std::endl;
         

         t = t0();
         u.resize(u0().size());
         k1.resize(u0().size());
         k2.resize(u0().size());
         k3.resize(u0().size());
         k4.resize(u0().size());
         k5.resize(u0().size());

        // std::vector<Type> Dt(u0().size(), dt());
         
         for(auto i=0 ; i < u0().size() ; i++)
            u[i] = u0()[i];                           //set Init Value  
         
         for(t = t0() ;  t <= tf() ; t += dt() )
         {
            std::cout << t << ' ' ;
            for(auto i=0 ; i< u.size() ; i++)                  // write Initial Value 
               std::cout << u[i] << " "; 
            std::cout << std::endl;

           for(auto j=0 ; j< u.size() ; j++)
           { 
            
             k1[j] = dt()*rhs.f[j](t                  , u                                          );   
             k2[j] = dt()*rhs.f[j](t+ 1.0/3.0*dt()    , u + k1[j]/3.0                              ); 
             k3[j] = dt()*rhs.f[j](t+ 1.0/3.0*dt()    , u + 1.0/6.0*(k1[j] + k2[j])                );
             k4[j] = dt()*rhs.f[j](t+ 1.0/2.0*dt()    , u + 1.0/8.0*(k1[j] + 3.0*k3[j] )           ); 
             k5[j] = dt()*rhs.f[j](t+dt()             , u + 1.0/2.0*(k1[j] - 3.0*k3[j] + 4.0*k4[j])); 
                  

            // 5-th order Runge Kutta (Fehlberg ) method 
              u[j] += 1.0/6.0 * ( k1[j] + 4.0*k4[j] + 1.0*k5[j] ) ;  
            }
            
         } 
         std::cout << "... Done " << std::endl;  
    

}

  }//ode
 }//numeric
}//mg 
# endif
