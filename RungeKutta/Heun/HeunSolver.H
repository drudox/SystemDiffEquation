# ifndef __HEUN_SOLVER_H__
# define __HEUN_SOLVER_H__

# include "../RungeKutta.H"
# include "../../rhsODEproblem.H"

namespace mg {
                namespace numeric {
                                     namespace odesystem {


/*-------------------------------------------------------------------|  
 *                                                                   |   
 *   Solve the differential problem :  y' = f(t,y)                   |
 *                                                 _______________   | 
 *   using Runge-Kutta Scheme (2 order accuracy)  >> Heun Method <<  |
 *                                                 ---------------   |
 *   @Marco Ghiani October 2017  -  Glasgow UK                       |
 *                                                                   |
 -------------------------------------------------------------------*/



template <typename Type>
class HeunSolver 
                  : public RungeKutta<Type>
{
  

  public:  

      HeunSolver(const rhsODEProblem<Type> & that) noexcept  :
                                                                RungeKutta<Type>{that} 
                  {}
      
      virtual ~HeunSolver() = default ;
      
      using OdeSolver<Type>::rhs;
      
      void solve(const std::string filename) override final  ;
      void solve() noexcept override final                   ;
     
    
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



template <typename Type>
inline void HeunSolver<Type>::solve(const std::string filename) {

      std::ofstream f(filename, std::ios::out );
      
      if(!f)
      {     
         std::cerr << "Error opening file " << filename << std::endl ;  
      }
      else
      {
         
         //t.at(0) = t0();
         
          u.resize(u0().size());
         k1.resize(u0().size());
         k2.resize(u0().size());
         
         for(auto i=0 ; i< u.size(); i++)
             u[i] = u0()[i];              // set I.V. 
            
         std::cout << "Running Heun (RK -2nd ord) Solver" << std::endl;
         
         for( t=t0() ; t <= tf() ; t+= dt() )
         {
            f << t << " " ;
            for(auto i=0 ; i< u.size(); i++)
               f << u[i] << " " ;
            f << std::endl;    

            for(auto i=0 ; i< u.size(); i++)
            {
               k1 = rhs.f[i](t       , u             );   
               k2 = rhs.f[i](t +dt() , u + dt()*k1[i]); 
            
               u[i] += dt()/2 *(k1[i]+k2[i]) ;  
            }
         } 
         std::cout << "... Done " << std::endl;  
      }
      f.close();
}

template <typename Type>
inline void HeunSolver<Type>::solve() noexcept
{
         
          u.resize(u0().size());
         k1.resize(u0().size());
         k2.resize(u0().size());
         
         for(auto i=0 ; i< u.size() ; i++)
             u[i] = u0()[i];              // set I.V. 
            
         std::cout << "Running Heun (RK -2nd ord) Solver" << std::endl;
         
         for( t=t0() ; t < tf() ; t+= dt() )
         {
            std::cout << t << " " ;
            for(auto i=0 ; i< u.size(); i++)
               std::cout << u[i] << " " ;
            std::cout << std::endl;    

            for(auto i=0 ; i< u.size(); i++)
            {
               k1 = rhs.f[i]( t       , u           );   
               k2 = rhs.f[i]( t +dt() , u+dt()*k1[i]); 
            
               u[i] += dt()/2 *(k1[i]+k2[i]) ;  
            }
         } 
         std::cout << "... Done " << std::endl;  
   
}
  
  }//ode
 }//numeric
}//mg 
# endif
