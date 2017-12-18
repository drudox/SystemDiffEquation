# ifndef __LEAPFROG_SOLVER_H__
# define __LEAPFROG_SOLVER_H__

# include "MultiStep.H"
# include "../rhsODEproblem.H"

namespace mg {
                namespace numeric {
                                     namespace odesystem {


/*-----------------------------------------------------------------------------
 *    Compute semi-Implicit solution  of a given (ODE) RHS problem using 
 *     ____________ 
 *    >> LEAP-FROG << scheme ( explicit 2 step - midpoint )
 *    -------------
 *    @Marco Ghiani Dec 2017, Glasgow UK
 *
 -----------------------------------------------------------------------------*/




template<typename Type = double>
class LeapFrogSolver :                                      // LEAP-FROG SOLVER
                                 public   MultiStep<Type> 
{
      
    public:  
      LeapFrogSolver(const rhsODEProblem<Type> & that) noexcept :
                                                                    MultiStep<Type>{that} 
                  {}
      
      virtual ~LeapFrogSolver() = default ;

      using OdeSolver<Type>::rhs;

      void solve(const std::string filename) override final ;
      void solve() noexcept override final                  ;
//
//
  private:

      using OdeSolver<Type>::t   ;
      using OdeSolver<Type>::u   ;
      
      using MultiStep<Type>::u_p1 ;
      using MultiStep<Type>::u_m1 ;
      
      using MultiStep<Type>::k1  ;
      using MultiStep<Type>::k2  ;
      
      using OdeSolver<Type>::dt  ; 
      using OdeSolver<Type>::t0  ;
      using OdeSolver<Type>::tf  ;
      using OdeSolver<Type>::u0  ;
      
      using OdeSolver<Type>::Ns  ;
      
      

};

//------------------  Implementation (to be put into .cpp file)   -----------------  //


template<typename Type>
inline void LeapFrogSolver<Type>::solve(const std::string filename) {
      
      std::ofstream f(filename, std::ios::out );
      
      if(!f)
      {     
         throw std::runtime_error(">> Error opening File in FWD euler solver <<" );
      }
      else
      {
      
         
         u.resize(u0().size());
         u_p1.resize(u0().size());
         u_m1.resize(u0().size());
         k1.resize(u0().size());
         k2.resize(u0().size());
 
         f << t0() << " " ;
         for(auto i=0 ; i< u.size(); i++ )   
         {
            u_m1[i] = u0()[i] ;                    // set initial value 
            f << u_m1[i] << " " ;                   // and write its                         
         }
         f << std::endl;
           
         //    
         //std::valarray<Type>  k1(u.size());
         //std::valarray<Type>  k2(u.size());
         
         for(auto i=0 ; i < u.size() ; i++)
         {
            k1[i] = rhs.f[i]( t          , u_m1 );
            k2[i] = rhs.f[i]( t + dt()/2 , u_m1 + k1[i]*dt()/2 );
            
                                             // initiation first point with 
            u[i] =  u_m1[i] + dt() * k2[i] ;  // Runge Kutta 2nd order PREDICTOR
         }
         
                                             
         
         std::cout << "Running LeapFrog (Leap-Frog) Solver" << std::endl;

         for(t=t0()+dt() ; t <= tf() ; t+=dt() )
         {
            f << t << " " ;
            for(auto i=0; i < u.size() ; i++)
              f << u[i] << " " ;
            f << std::endl;
            
            for(auto i=0 ; i < u.size() ; i++)
            {
               u_p1[i] = u_m1[i] + 2*dt() * rhs.f[i]( t , u );  // leap-frog 
            
            
               u_m1[i] = u[i];
                 u[i] = u_p1[i];
            }

         } 
         std::cout << "... Done " << std::endl;  
         f.close();
      } 
}


template<typename Type>
inline void LeapFrogSolver<Type>::solve() noexcept 
{

   
         u.resize(u0().size());
         u_p1.resize(u0().size());
         u_m1.resize(u0().size());
         
         std::cout << t0() << " " ;
         for(auto i=0 ; i< u.size() ; i++ )   
         {
            u_m1[i] = u0()[i] ;                    // set initial value 
            std::cout << u_m1[i] << " " ;                   // and write its                         
         }
         std::cout << std::endl;
           
         //    
         
         for(auto i=0 ; i < u.size() ; i++)
         {
            k1[i] = rhs.f[i]( t          , u );
            k2[i] = rhs.f[i]( t + dt()/2 , u + k1[i]*dt()/2 );
            
                                             // initiation first point with 
            u[i] =  u_m1[i] + dt() * k2[i] ;  // Runge Kutta 2nd order PREDICTOR
         }
         
                                             
         
         std::cout << "Running LeapFrog (Leap-Frog) Solver" << std::endl;

         for(t=t0()+dt() ; t < tf() ; t+=dt() )
         {
            std::cout << t << " " ;
            for(auto i=0; i < u.size() ; i++)
              std::cout << u[i] << " " ;
            std::cout << std::endl;
            
            for(auto i=0 ; i < u.size() ; i++)
            {
               u_p1[i] = u_m1[i] + 2*dt() * rhs.f[i](t , u );  // leap-frog 
            
               u_m1[i] = u[i];
                 u[i]  = u_p1[i];
            }

         } 
         std::cout << "... Done " << std::endl;  
       



}
  
  }//ode
 }//numeric
}//mg
# endif
