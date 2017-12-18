# ifndef __RUNGEKUTTA4_SOLVER_H__
# define __RUNGEKUTTA4_SOLVER_H__

//# include "OdeSolver.H"
# include "rhsODEproblem.H"
# include "RungeKutta.H" 

namespace mg { 
                namespace numeric {
                                    namespace odesystem {


/*-------------------------------------------------------------------------------
 *    
 *    Perform Runge-Kutta (4th order accuracy) solution of (ODE) RHS problem
 *    dy/dt = f(y,t)
 *
 *    @Marco Ghiani Dec 2017, Glasgow UK
 *
 ------------------------------------------------------------------------------*/



template<typename Type= double>
class RungeKutta4Solver 
                         :   public  RungeKutta<Type>    
{
      
    public:  
      RungeKutta4Solver(const rhsODEProblem<Type> & that) noexcept :
                                                                        RungeKutta<Type>{that} 
                  {}
      
      virtual ~RungeKutta4Solver() = default;
      

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
 

};


template<typename Type>
inline void RungeKutta4Solver<Type>::solve(const std::string filename)  {
      
      std::ofstream f(filename, std::ios::out );
      
      if(!f)
      {     
         std::string mess = "Error opening file " + filename + " in RungeKutta4-Solver " ;  
         throw std::runtime_error(mess.c_str());   
      }
      else
      {
         std::cout << "Running Runge-Kutta 4th order Solver" << std::endl;
         

         t = t0();
          u.resize(u0().size());
         k1.resize(u0().size());
         k2.resize(u0().size());
         k3.resize(u0().size());
         k4.resize(u0().size());

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
            
                k1[j] = rhs.f[j](t        , u                  );   
                k2[j] = rhs.f[j](t+dt()/2 , u + dt()/2* k1[j]  ); 
                k3[j] = rhs.f[j](t+dt()/2 , u + dt()/2* k2[j]  );
                k4[j] = rhs.f[j](t+dt()   , u + dt() *  k3[j]  ); 
            
                u[j] += dt()/6 *(k1[j]+ 2*k2[j] + 2*k3[j] +k4[j]) ;  
            }
            
         } 
         std::cout << "... Done " << std::endl;  
      
      f.close();
      }
}


template<typename Type>
inline void RungeKutta4Solver<Type>::solve() noexcept 
{
      std::cout << "Running Runge-Kutta 4th order Solver" << std::endl;
         

      t = t0();
      u.resize(u0().size());
      k1.resize(u0().size());
      k2.resize(u0().size());
      k3.resize(u0().size());
      k4.resize(u0().size());
      
      std::vector<Type> Dt(u0().size(), dt());


      for(auto i=0 ; i < u0().size() ; i++)
          u[i] = u0()[i];                           //set Init Value  
         
      for(t = t0() ;  t <= tf() ; t += dt() )
      {
         std::cout << t << ' ' ;
         for(auto i=0 ; i< u.size() ; i++)                  // write Initial Value 
               std::cout << u[i] << std::endl; 
         std::cout << std::endl;

           for(auto j=0 ; j< u.size(); j++)
           { 
            
                k1[j] = rhs.f[j](t        , u                );   
                k2[j] = rhs.f[j](t+dt()/2 , u + dt()/2*k1[j] ); 
                k3[j] = rhs.f[j](t+dt()/2 , u + dt()/2*k2[j] );
                k4[j] = rhs.f[j](t+dt()   , u + dt()*k3[j]   ); 
            
                u[j] += dt()/6 *(k1[j]+ 2*k2[j] + 2*k3[j] +k4[j]) ;  
            }
            
       } 
       std::cout << "... Done " << std::endl;  
      
 

}

  }//ode
 }//numeric
}//mg 
# endif
