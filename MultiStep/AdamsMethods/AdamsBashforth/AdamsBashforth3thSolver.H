# ifndef __ADAMS_BASHFORTH_3TH_SOLVER_H__
# define __ADAMS_BASHFORTH_3TH_SOLVER_H__

# include "../AdamsMethods.H"
# include "../../../rhsODEproblem.H"

namespace mg {
                namespace numeric {
                                     namespace odesystem {


/**-----------------------------------------------------------------------------
 *
 *\ @brief \class AdamsBashforth3thSolver
 *\ Perform linear multi-step ( 3 ) solution of a given (ODE) RHS problem 
 * 
 *
 *\ @Marco Ghiani Dec 2017, Glasgow UK
 *
 *
 -----------------------------------------------------------------------------*/




template<typename Type = double>
class AdamsBashforth3thSolver :
                                 public AdamsMethods<Type>
{
      
    public:  

      AdamsBashforth3thSolver(const rhsODEProblem<Type> & that) noexcept :
                                                                            AdamsMethods<Type>{that} 
                  {}
      
      virtual ~AdamsBashforth3thSolver() = default ;

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
      using AdamsMethods<Type>::up1 ;

      using AdamsMethods<Type>::K1 ;
      using AdamsMethods<Type>::K2 ;
      
 

};

//------------------  Implementation (to be put into .cpp file)   -----------------  //


template<typename Type>
inline void AdamsBashforth3thSolver<Type>::solve(const std::string filename) {
      
      std::ofstream f(filename, std::ios::out );
      
      if(!f)
      {     
         throw std::runtime_error(">> Error opening File in ADAMS-BASHFORTH  solver <<" );
      }
      else
      {
         std::cout << "Running Adams-Bashforth 3th order (3 step) Solver" << std::endl;
      
         u.resize(u0().size());
         um2.resize(u0().size());
         um1.resize(u0().size());
         up1.resize(u0().size());
         
         K1.resize(u0().size());
         K2.resize(u0().size());


         t = t0();

         f << t0() << ' ' ;
         for(auto i=0 ; i < u.size() ; i++)
         {
            um2[i] = u0()[i] ;                  // set initial value 
            f << um2[i] << ' ' ;                // and write there to file
         }
         f << std::endl;
         
      
         // using for the trigger points  (start up the method)
         //
          
         for(auto i=0 ; i< u.size() ; i++)
         {  
            K1[i] = rhs.f[i]( t        ,  um2 );
            K2[i] = rhs.f[i]( t + dt() ,  um2 + dt()*K1[i] );
            // compute first point 
            um1[i] =  um2[i] + dt()/2 *(K1[i] + K2[i]) ;  // Rk 2nd order 
         }
         
         f << t+dt() << " " ;
         for(auto i=0 ; i < u.size() ; i++)
            f << um1[i] << "  " ;
         f << std::endl;   


         for(auto i=0; i < u.size() ; i++)
         {
            K1[i] = rhs.f[i]( t +  dt() , um1 );
            K2[i] = rhs.f[i]( t + 2*dt(), um1+ dt()*K1[i] );
            // compute second point
            u[i]   = um1[i] + dt()/2 *(K1[i] + K2[i]);  // rk 2nd order 
         }

         for(t=t0()+dt()+dt() ; t <= tf() ; t += dt() )
         {
            f << t << "  " ;
            for(auto i=0 ; i < u.size() ; i++)
                  f << u[i] << " " ;
            f << std::endl; 
            
            for(auto i=0 ; i< u.size(); i++)
            { 
               up1[i] = u[i] + dt()/12 *(23 * rhs.f[i](t , u) - 16* rhs.f[i](t-dt() , um1) + 
                                           5* rhs.f[i](t-dt()-dt() , um2)  )  ;
            }
            for(auto i=0 ; i < u.size() ; i++)
            {
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
inline void AdamsBashforth3thSolver<Type>::solve() noexcept 
{
       
         std::cout << "Running Adams-Bashforth (3 step) Solver" << std::endl;
      
         u.resize(u0().size());
         um2.resize(u0().size());
         um1.resize(u0().size());
         up1.resize(u0().size());
       
         // k** using for the trigger points  (start up the method)
         K1.resize(u.size());
         K2.resize(u.size());
        
 
   
         t = t0();

         std::cout << t0() << ' ' ;
         for(auto i=0 ; i < u.size() ; i++)
         {
            um2[i] = u0()[i] ;                  // set initial value 
            std::cout << um2[i] << ' ' ;        // and print on cout 
         }
         std::cout << std::endl;
         
      
         for(auto i=0 ; i< u.size() ; i++)
         {  
            K1[i] = rhs.f[i]( t        ,  um2 );
            K2[i] = rhs.f[i]( t + dt() ,  um2 + dt()*K1[i] );
            // compute first point 
            um1[i] =  um2[i] + dt()/2 *(K1[i] + K2[i]) ;  // Rk 2nd order 
         }
         
         std::cout << t+t0() << " " ;
         for(auto i=0 ; i < u.size() ; i++)
            std::cout << um1[i] << "  " ;
         std::cout << std::endl;   


         for(auto i=0; i < u.size() ; i++)
         {
            K1[i] = rhs.f[i]( t +  dt() , um1 );
            K2[i] = rhs.f[i]( t + 2*dt(), um1+ dt()*K1[i] );
            // compute second point
            u[i]   = um1[i] + dt()/2 *(K1[i] + K2[i]);  // rk 2nd order 
         }

         for(t=t0()+dt()+dt() ; t <= tf() ; t += dt() )
         {
            std::cout << t << "  " ;
            for(auto i=0 ; i < u.size() ; i++)
                 std::cout << u[i] << " " ;
            std::cout << std::endl; 
            
            for(auto i=0 ; i< u.size(); i++)
            { 
               up1[i] = u[i] + dt()/12 *(23 * rhs.f[i](t , u) - 16* rhs.f[i](t-dt() , um1) + 
                                           5* rhs.f[i](t-dt()-dt() , um2)  )  ;
            }
            for(auto i=0 ; i < u.size() ; i++)
            {
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
