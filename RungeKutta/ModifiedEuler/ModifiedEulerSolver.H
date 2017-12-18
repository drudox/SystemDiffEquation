#ifndef __MODIFIED_EULER_H__
#define __MODIFIED_EULER_H__


# include "rhsODEproblem.H"
# include "RungeKutta.H"

namespace mg { 
               namespace numeric {
                                    namespace odesystem {

/*-----------------------------------------------------------------------
 *    compute the solution of RHS (ODE) problem  y' = f(t,y) 
 *    using (Runge Kutta 2th order accuracy) 
 *    
 *    Modified - explicit Euler       
 *    
 *    @Marco Ghiani Dec. 2017 Glasgow 
 *
 ----------------------------------------------------------------------*/




template <typename Type>
class ModifiedEulerSolver 
                          : public RungeKutta<Type>
{


   public:   

    ModifiedEulerSolver(const rhsODEProblem<Type>& that) noexcept :
                                                                    RungeKutta<Type>{that} 
                    {}
    
    virtual ~ModifiedEulerSolver() = default;
    

    using OdeSolver<Type>::rhs;


    void solve(const std::string filename ) override final;
    void solve() noexcept override final                  ;

//--
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
      
};


//----------------- Implementation (to be put into .cpp file) ----------

template <typename Type>
inline void ModifiedEulerSolver<Type>::solve(const std::string fname)  
{
   std::ofstream f(fname , std::ios::out);

   if(!f)
   {
      std::string mess = "Error opening file " + fname + "in Modified Euler Solver " ; 
      throw std::runtime_error(mess.c_str());
   }
   
   //t = t0();
    u.resize(u0().size());   
   k1.resize(u0().size());   
   k2.resize(u0().size());   
   
   for(auto i=0 ; i < u.size() ; i++)
       u[i] = u0()[i] ;                   // setting initial value
      


      
   std::cout << "Running Modified Euler (RK -2nd ord) Solver" << std::endl;
   for(t = t0() ; t <= tf() ; t+=dt()  )
   {   
       f << t << ' ' ;
       for(auto i=0; i< u.size(); i++)
          f  << u[i] << "  ";
       f << std::endl;

       for(auto i=0; i < u.size() ; i++)
       {
          k1[i] = rhs.f[i](t , u );
          k2[i] = rhs.f[i](t + dt()/2 , u + k1[i]*dt()/2 );
      
          u[i] += dt() * k2[i] ;
       }  
   }
   f.close();
   std::cout << "... Done" << std::endl;
}


template <typename Type>
inline void ModifiedEulerSolver<Type>::solve() noexcept 
{
  //t = t0();
    u.resize(u0().size());   
   k1.resize(u0().size());   
   k2.resize(u0().size());   
   
   for(auto i=0 ; i < u.size() ; i++)
       u[i] = u0()[i] ;                   // setting initial value
      


      
   std::cout << "Running Modified Euler (RK -2nd ord) Solver" << std::endl;
   for(t = t0() ; t < tf() ; t+=dt()  )
   {   
       std::cout << t << ' ' ;
       for(auto i=0; i< u.size(); i++)
            std::cout << u[i] << std::endl;
       std::cout << std::endl;

       for(auto i=0; i < u.size() ; i++)
       {
          k1[i] = rhs.f[i](t , u );
          k2[i] = rhs.f[i](t + dt()/2 , u + k1[i]*dt()/2 );
      
          u[i] += dt() * k2[i] ;
       }  
   }
   std::cout << "... Done" << std::endl;
}
  
  }//ode
 }//numeric
}//mg 
# endif
