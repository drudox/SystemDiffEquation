# ifndef BACKWARD_EULER_SOLVER_H
# define BACKWARD_EULER_SOLVER_H

# include "Euler.H"
# include "../rhsODEproblem.H"

namespace mg { 
               namespace numeric {
                                    namespace odesystem {


/*----------------------------------------------------------------------*
 *    
 *    Compute Implicit Euler , solved the non-linear algebric equations
 *    using Newton-Rapson 
 *
 *    Solution of a given (ODE) RHS problem :
 *    y' = f(t,y); 
 *
 *    @Marco Ghiani , Dec 2017 , Glasgow UK
 *
 -----------------------------------------------------------------------*/



template <typename Type = double>
class BackwardEulerSolver : 
                              public Euler<Type>
{
      
   public:  
      BackwardEulerSolver(const rhsODEProblem<Type> & that) :
                                                               Euler<Type>{that} 
                  {}

      virtual ~BackwardEulerSolver() = default ;
      
      using OdeSolver<Type>::rhs;

      void solve(std::string filename) override final ;
      void solve()  noexcept override final ;
   
   private:
      
      using OdeSolver<Type>::t  ;
      using OdeSolver<Type>::u  ;
      using OdeSolver<Type>::dt ; 
      using OdeSolver<Type>::t0 ;
      using OdeSolver<Type>::tf ;
      using OdeSolver<Type>::u0 ;
      
      using OdeSolver<Type>::Ns ;
      
      using OdeSolver<Type>::toll ;

      std::valarray<Type> uOld ;
      std::valarray<Type> uNew ;
      std::valarray<Type> err  ; 

      
};

//------------------  Implementation (to be put into .cpp file) -------------------- //


template <typename Type>
inline void BackwardEulerSolver<Type>::solve(std::string filename)  {
      
      std::ofstream f(filename, std::ios::out );
      
      if(!f)
      {     
         std::string mess = "Error opening file " + filename + "in BackwardEulerSolver" ;  
         throw std::runtime_error(mess.c_str());         
      }
      else
      {
         std::cout << "Running BackwardEuler Solver" << std::endl;
         
         
         t = t0();
         u.resize(u0().size());

         uOld.resize(u0().size());
         uNew.resize(u0().size());
         err.resize(u0().size());
         
         for(auto i=0; i < u0().size() ; i++ )        // setting initail value
            u[i] = u0()[i] ;
  
         f << t << ' ' ;
         for(auto i=0 ; i < u0().size() ; i++ )       // write initial value 
            f << u[i] << std::endl; 
         f << std::endl;
            
         auto max_it = 100;
         for(t=t0()+dt() ; t <= tf() ; t+=dt() )
         {
            
            for(auto j=0 ; j< u.size() ; j++)
            {   
              
              uOld[j] = u[j] + dt() *rhs.f[j](t+dt(), u) ;  // semi implicit-predictor 
              
              err[j] = 1.0;

              auto k=1;
              while(err[j] > toll ) //&& j < max_it)
              {
                  
                //  uNew[j] = u[j] + dt()*rhs.f[j](t,uOld);               
                  
                  uNew[j] = uOld[j] - ( uOld[j] - (u[j] + dt() * rhs.f[j](t,uOld) ) ) /
                                 ( 1- dt() * rhs.dfdt(j, t , uOld) ) ;
                                                                               
                  err[j] = fabs(uNew[j] - uOld[j] );
               
                  uOld[j] = uNew[j];
               }
               k++;
            }      
            
            for(auto j=0 ; j < u.size() ; j++)
               u[j] = uNew[j] ;
            
            f << t << ' ' ;
            for(auto j=0 ; j< u.size() ; j++ ){
                   f << u[j] <<  " " ;
            }       
            f <<std::endl ;
         } 
         std::cout << "... Done " << std::endl;  
      f.close();
      }
}


template <typename Type>
inline void BackwardEulerSolver<Type>::solve() noexcept 
{
         std::cout << "Running BackwardEuler Solver" << std::endl;
         
         t = t0();

         u.resize(u0().size());
         uOld.resize(u0().size());
         err.resize(u0().size());
         
         
         for(auto i=0; i < u0().size() ; i++ )
            u[i] = u0()[i] ;
  
         std::cout << t << ' ' ;
         for(auto i=0 ; i < u0().size() ; i++ )
            std::cout << u[i] << std::endl; 
         std::cout << std::endl;
            
         auto max_it = 100;
         for(t=t0()+dt() ; t < tf() ; t+=dt())
         {
            
            for(auto j=0 ; j< u.size() ; j++)
            {   
              
      //        uOld[j] = u[j] + dt() *rhs.f[j](t, u) ;  
              err[j] = 1.0;
            
              auto k=1;
              while(err[j] > toll ) //&& k < max_it)
              {   

        //          uNew[j] = u[j] + dt()*rhs.f[j](t,uOld);               
                  uNew[j] = uOld[j] - ( uOld[j] - (u[j] + dt() * rhs.f[j](t,uOld) ) ) /
                                 ( 1- dt() * rhs.dfdt(j, t , uOld) ) ;
                
                  err[j] = fabs(uNew[j] - uOld[j] );
               
                  uOld[j] = uNew[j];
               k++;
               }
            }      
            
            for(auto j=0 ; j < u.size() ; j++)
               u[j] = uNew[j] ;
            
            std::cout << t << ' ' ;
            for(auto j=0 ; j< u.size() ; j++ ){
                   std::cout << u[j] <<  " " ;
            }       
            std::cout <<std::endl ;
         } 
         std::cout << "... Done " << std::endl;  
}

  }//ode
 }//numeric
}//mg 
# endif
