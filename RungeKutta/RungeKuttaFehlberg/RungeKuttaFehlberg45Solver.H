# ifndef __RKF45_SOLVER_H__
# define __RKF45_SOLVER_H__



/*
 *
 *
 */


# include "../RungeKutta.H"
# include "../../rhsODEproblem.H"

namespace mg { 
               namespace numeric {
                                    namespace odesystem {

/**------------------------------------------------------------------
 * \class RungeKuttaFehlberg45Solver  
 * \brief Adaptative step size, checking the step error (every step) 
 *
 * \author  Marco Ghiani
 * \date       Jan 2018 
 * \place      Glasgow
 *
 -------------------------------------------------------------------*/


template <typename Type = double> 
class RungeKuttaFehlberg45Solver      
                                  : public RungeKutta<Type>
{

    public:

 //-- std constructor
      RungeKuttaFehlberg45Solver (const rhsODEProblem<Type> & that ) noexcept : 
                                                            RungeKutta<Type>{that} , 
                                                            h0{that.dt()}      
                  {}

      virtual ~RungeKuttaFehlberg45Solver() = default ;
      
      using OdeSolver<Type>::rhs ;


      void solve(const std::string fname) override final ;

      void solve() noexcept override final ;
     

//--     
    private:   

      using OdeSolver<Type>::t ;
      using OdeSolver<Type>::u  ;
      using OdeSolver<Type>::dt ; 
      using OdeSolver<Type>::t0 ;
      using OdeSolver<Type>::tf ;
      using OdeSolver<Type>::u0 ;
      

      using OdeSolver<Type>::Ns ;
      
      using RungeKutta<Type>::u1 ;
      using RungeKutta<Type>::u2 ;

      //Type u2 ;
      
      using OdeSolver<Type>::toll ;
      using OdeSolver<Type>::setTimeStep ;

      using RungeKutta<Type>::k1 ;
      using RungeKutta<Type>::k2 ;
      using RungeKutta<Type>::k3 ;
      using RungeKutta<Type>::k4 ;
      using RungeKutta<Type>::k5 ;
      using RungeKutta<Type>::k6 ;
      
      using RungeKutta<Type>::s;
      using RungeKutta<Type>::h;
      using RungeKutta<Type>::delta;
      using RungeKutta<Type>::R;
      


      const Type h0 ;
     // Type delta ;
      // Type R     ;
      using RungeKutta<Type>::stepToll ;// = 1.0e-6 ;
      using RungeKutta<Type>::epsilon  ;//  = 1.0e-6 ;

};

template<typename Type> 
inline void RungeKuttaFehlberg45Solver<Type>::solve(std::string filename) {
      
   std::ofstream f(filename , std::ios::out);

   if(!f)
   {
      std::string mess = "Error opening file " + filename + " in RKF45 Solver " ; 
      throw std::runtime_error(mess.c_str());      
   }
   else
   {

      std::cout << "Running Runge-Kutta Fehlberg 4-5 order Solver " << std::endl; 
      
      t = t0();
      u.resize(u0().size());
      u1.resize(u0().size());
      u2.resize(u0().size());
      
      s.resize(u0().size());
      //h.resize(u0().size());
      delta.resize(u0().size());
      R.resize(u0().size());
      
      k1.resize(u0().size());
      k2.resize(u0().size());
      k3.resize(u0().size());
      k4.resize(u0().size());
      k5.resize(u0().size());
      k6.resize(u0().size());

      for (auto i=0 ; i< u0().size() ; i++ )
      {
         u[i]  = u0()[i] ;
         u1[i] = u0()[i] ;
         u2[i] = u0()[i] ;

      }

      h = dt() ;      

      for(t=t0() ; t<= tf() ;  )
      {
            f << t << ' ' ;
            for(auto i=0 ; i< u.size() ; i++)                  // write Initial Value 
               f << u[i] << " "; 
            f << std::endl;
           

           for(auto j=0 ; j< u.size() ; j++)
           { 
            
             k1[j] = h*rhs.f[j](t                  , u                                             );   
             k2[j] = h*rhs.f[j](t+dt()/4           , u + k1[j]/4.0                                 ); 
             k3[j] = h*rhs.f[j](t+ 3.0/8.0*dt()    , u + 3.0/32.0*k1[j] + 9.0/32.0* k2[j]          );
             k4[j] = h*rhs.f[j](t+ 12.0/13.0* dt() , u + 1932.0/2197.0* k1[j] -  7200.0/2197.0*k2[j] 
                                                       + 7296.0/2197.0* k3[j]                      ); 
             k5[j] = h*rhs.f[j](t+dt()             , u + 439.0/216.0  * k1[j] - 8.0*k2[j]  
                                                       + 3680.0/513.0 * k3[j] - 845.0/4104.0*k4[j] ); 
             k6[j] = h*rhs.f[j](t+dt()/2.0         , u - 8.0/27.0 * k1[j] + 2.0*k2[j] 
                                                          - 3544.0/2565.0* k3[j] + 1859.0/4104.0*k4[j]
                                                                                   - 11.0/40.0    *k5[j]); 
                  

                        // 4-th order Runge Kutta (Fehlberg) method
             u1[j] = u[j] + (25.0/216.0*k1[j] + 1408.0/2565.0*k3[j] + 2197.0/4101.0*k4[j]-1.0/5.0*k5[j]);   
                       
                       // 5-th order Runge Kutta (Fehlberg ) method 
             u2[j] = u[j] + (16.0/135.0*k1[j] + 6656.0/12825.0 *k3[j] + 28561.0/56430.0* k4[j] 
                                                 - 9.0/50.0* k5[j] + 2.0/55.0* k6[j]) ;  
            

             R[j] = fabs(u1[j]-u2[j])/ h ; //relative Error 
            
             s[j] = 0.8409*(pow(epsilon/R[j],(1./4)));  // adaptative Time-step
                 
                 
             //u[j] = u_5[j] ;
             
             if(R[j] <= epsilon) 
             {
                 t = t+h;
                 u[j] = u1[j] ;
                 h = s.min() * h ; 
             }
             else
             {
                h = s.min() *h ;   
             }



            }
            

             

      }    
      std::cout << "... done" << std::endl;
   }   
}

template<typename Type> 
inline void RungeKuttaFehlberg45Solver<Type>::solve() noexcept {
      
 //  std::ofstream f(filename , std::ios::out);

      std::cout << "Running Runge-Kutta Fehlberg 4-5 order Solver " << std::endl; 
      
      t = t0();
      u.resize(u0().size());
      u1.resize(u0().size());
      u2.resize(u0().size());
      
      s.resize(u0().size());
      //h.resize(u0().size());
      delta.resize(u0().size());
      R.resize(u0().size());
      
      k1.resize(u0().size());
      k2.resize(u0().size());
      k3.resize(u0().size());
      k4.resize(u0().size());
      k5.resize(u0().size());
      k6.resize(u0().size());

      for (auto i=0 ; i< u0().size() ; i++ )
      {
         u[i]  = u0()[i] ;
         u1[i] = u0()[i] ;
         u2[i] = u0()[i] ;

      }

      h = dt() ;      

      for(t=t0() ; t<= tf() ;  )
      {
            std::cout << t << ' ' ;
            for(auto i=0 ; i< u.size() ; i++)                  // write Initial Value 
               std::cout << u[i] << " "; 
            std::cout << std::endl;
           

           for(auto j=0 ; j< u.size() ; j++)
           { 
            
             k1[j] = h*rhs.f[j](t                  , u                                             );   
             k2[j] = h*rhs.f[j](t+dt()/4           , u + k1[j]/4.0                                 ); 
             k3[j] = h*rhs.f[j](t+ 3.0/8.0*dt()    , u + 3.0/32.0*k1[j] + 9.0/32.0* k2[j]          );
             k4[j] = h*rhs.f[j](t+ 12.0/13.0* dt() , u + 1932.0/2197.0* k1[j] -  7200.0/2197.0*k2[j] 
                                                       + 7296.0/2197.0* k3[j]                      ); 
             k5[j] = h*rhs.f[j](t+dt()             , u + 439.0/216.0  * k1[j] - 8.0*k2[j]  
                                                       + 3680.0/513.0 * k3[j] - 845.0/4104.0*k4[j] ); 
             k6[j] = h*rhs.f[j](t+dt()/2.0         , u - 8.0/27.0 * k1[j] + 2.0*k2[j] 
                                                          - 3544.0/2565.0* k3[j] + 1859.0/4104.0*k4[j]
                                                                                   - 11.0/40.0    *k5[j]); 
                  

                        // 4-th order Runge Kutta (Fehlberg) method
             u1[j] = u[j] + (25.0/216.0*k1[j] + 1408.0/2565.0*k3[j] + 2197.0/4104.0*k4[j]-1.0/5.0*k5[j]);   
                       
                       // 5-th order Runge Kutta (Fehlberg ) method 
             u2[j] = u[j] + (16.0/135.0*k1[j] + 6656.0/12825.0 *k3[j] + 28561.0/56430.0* k4[j] 
                                                 - 9.0/50.0* k5[j] + 2.0/55.0* k6[j]) ;  
            

             R[j] = fabs(u1[j]-u2[j])/ h ; //relative Error 
            
             s[j] = 0.8409*(pow(epsilon/R[j],(1./4)));  // adaptative Time-step
                 
                 
             //u[j] = u_5[j] ;
             
             if(R[j] <= epsilon) 
             {
                 t = t+h;
                 u[j] = u1[j] ;
                 h = s.min() * h ; 
             }
             else
             {
                h = s.min() *h ;   
             }



            }
            

             

          

   }   
}


      

  }//odesystem
 }// numeric
}//mg 
# endif 
