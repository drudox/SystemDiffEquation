# ifndef __ADAMS_MOULTON_4ORD_3STEP_SOLVER_H__
# define __ADAMS_MOULTON_4ORD_3STEP_SOLVER_H__


# include "../AdamsMethods.H"
# include "../../../rhsODEproblem.H"


namespace mg {
                namespace numeric {
                                     namespace odesystem {


/**------------------------------------------------------------------------------------
 *
 * @class AdamsMoulton4thSolver
 * @brief class Adams Bashforth - Predictor 4th order - corrector 3th Solver
 * 
 *  Perform linear multi-step (2) scheme using PREDICTOR (Adams-bashforth 4 step) 
 *                                           & CORRECTOR (Adams-PC4th (imp) 3 step )  
 * 
 *
 * @author Marco Ghiani Dec 2017, Glasgow UK
 *
 *
 ------------------------------------------------------------------------------------*/



template<typename Type = double>
class AdamsMoulton4thSolver :
                                 public AdamsMethods<Type> 
{

//---      
    public:  

      AdamsMoulton4thSolver(const rhsODEProblem<Type> & that) noexcept :
                                                                          AdamsMethods<Type>{that} 
                  {}
      
      virtual ~AdamsMoulton4thSolver() = default ;

      using OdeSolver<Type>::rhs;

      void solve(const std::string filename) override final;
      void solve() noexcept override final                 ;
//
//--
  private:

      using OdeSolver<Type>::t  ;
      using OdeSolver<Type>::u  ;
      using OdeSolver<Type>::dt ; 
      using OdeSolver<Type>::t0 ;
      using OdeSolver<Type>::tf ;
      using OdeSolver<Type>::u0 ;
      
      using OdeSolver<Type>::toll ;
      
      using AdamsMethods<Type>::um1 ;
      using AdamsMethods<Type>::um2 ; 
      using AdamsMethods<Type>::um3 ; 
      using AdamsMethods<Type>::up1 ;
 //     using AdamsMethods<Type>::u   ;  
      
      
      using AdamsMethods<Type>::uPred;
      using AdamsMethods<Type>::uCorrOld;   // used during iteration
      using AdamsMethods<Type>::uCorr;
      
      using AdamsMethods<Type>::fPred;
      using AdamsMethods<Type>::fCorr;
      using AdamsMethods<Type>::fCorrOld;
      
      using AdamsMethods<Type>::K1;
      using AdamsMethods<Type>::K2;
      using AdamsMethods<Type>::K3;
      using AdamsMethods<Type>::K4;
      
      using AdamsMethods<Type>::error;


};

//------------------  Implementation (to be put into .cpp file)   -----------------  //


template<typename Type>
inline void AdamsMoulton4thSolver<Type>::solve(const std::string filename) {
      
    std::ofstream f(filename, std::ios::out );
      
    if(!f)
    {     
       throw std::runtime_error(">> Error opening File in ADAMS BASHFORTH-MOULTON  solver <<" );
    }
    else
    {
         std::cout << "Running Adams Bashforth (4 step) - Adams Moulton (3 step) Solver" << std::endl;
      
         u.resize(u0().size());
         um1.resize(u0().size());
         um2.resize(u0().size());
         um3.resize(u0().size());
         up1.resize(u0().size());
         
         uPred.resize(u0().size());
         uCorr.resize(u0().size());
         uCorrOld.resize(u0().size());
         
         fPred.resize(u0().size());
         fCorr.resize(u0().size());
         fCorrOld.resize(u0().size());    // using to iterate 
         
         K1.resize(u0().size());
         K2.resize(u0().size());
         K3.resize(u0().size());
         K4.resize(u0().size());
         
         error.resize(u0().size());

         t = t0();

         f << t << ' ' ;
         for(auto i=0 ; i < u.size() ; i++)
         {
            um3[i] = u0()[i] ;                  // set initial value 
            f << um3[i] << ' ' ;                // and write there to file
         }
         f << std::endl;
         
         
         
         for(auto i=0 ; i< u.size() ; i++)
         {  
            K1[i] = rhs.f[i]( t          ,  um3 );
            K2[i] = rhs.f[i]( t + dt()/2 ,  um3 + dt()/2 * K1[i] );
            K3[i] = rhs.f[i]( t + dt()/2 ,  um3 + dt()/2 * K2[i] );
            K4[i] = rhs.f[i]( t + dt()   ,  um3 + dt()   * K3[i] );
                                      // compute first point (start - up  Adams Bashforth )
            um2[i] =  um3[i] + dt()/6 *(K1[i] + 2*K2[i] + 2*K3[i] + K4[i]) ;  // Rk 4th order 
         }
            

         t+=dt();   
            
         f << t << ' ' ;
         for(auto i=0 ; i < u.size() ; i++)
            f << um2[i] << ' ' ;                // and write there to file
         f << std::endl;
    
         for(auto i=0 ; i< u.size() ; i++)
         {  
            K1[i] = rhs.f[i]( t          ,  um2 );
            K2[i] = rhs.f[i]( t + dt()/2 ,  um2 + dt()/2*K1[i] );
            K3[i] = rhs.f[i]( t + dt()/2 ,  um2 + dt()/2*K2[i] );
            K4[i] = rhs.f[i]( t + dt()   ,  um2 + dt()  *K3[i] );
            
                                      // compute second point (start - up  Adams Bashforth )
            um1[i] =  um2[i] + dt()/6 *(K1[i] + 2*K2[i] +2*K3[i] + K4[i] ) ;  // Rk 4th order 
         }

         t+=dt();
         
         f << t << ' ' ;
         for(auto i=0 ; i < u.size() ; i++)
            f << um1[i] << ' ' ;                // and write there to file
         f << std::endl;
         
         for(auto i=0; i < u.size() ; i++)
         {
            K1[i] = rhs.f[i]( t            , um1 );
            K2[i] = rhs.f[i]( t +   dt()/2 , um1+ dt()/2 *K1[i] );
            K3[i] = rhs.f[i]( t +   dt()/2 , um1+ dt()/2 *K2[i] );
            K4[i] = rhs.f[i]( t +   dt()   , um1+ dt()   *K3[i] );
            // compute second point
            u[i]   = um1[i] + dt()/6 *(K1[i] + 2*K2[i] + 2*K3[i] + K4[i]);  // rk 2nd order 
         }

     // MAIN LOOP 
         for(t = t+dt() ; t < tf() ; t += dt() )
         {
            f << t << "  " ;
            for(auto i=0 ; i < u.size() ; i++)
                  f << u[i] << " " ;
            f << std::endl; 


//------------ PREDICTOR      ADAMS-BASHFORTH 4TH ORDER (4 step)
//
            for(auto i=0 ; i< u.size(); i++)
            { 
                  uPred[i] = u[i] + dt()/24.0 *( 55.0 * rhs.f[i](t , u) 
                                               - 59.0 * rhs.f[i](t-dt() , um1) 
                                               + 37.0 * rhs.f[i](t-(2*dt()), um2 ) 
                                               -  9.0 * rhs.f[i](t-(3*dt()), um3 ) );
            

                  fPred[i] = rhs.f[i](t+dt() , uPred);     
            }

//------------ CORRECTOR     ADAMS-MOULTON 4TH ORDER (3 step)          
            
            for(auto i=0; i< u.size() ; i++){
                  error[i] = 1.0 ;        
                  fCorrOld[i] = fPred[i] ;
            }

            std::size_t iter = 0 ;
           
          for(auto i=0 ; i < u.size() ; i++)
          {
           while( error[i] >= toll )
           {
              uCorrOld[i] = u[i] + dt()/24.0*(   9.0* fCorrOld[i] 
                                              + 19.0* rhs.f[i](t,u) 
                                              -  5.0* rhs.f[i](t-dt(), um1 )
                                              +  1.0* rhs.f[i](t-(2*dt()), um2 )  );  
                  
              fCorrOld[i] = rhs.f[i](t+dt(), uCorrOld );

              uCorr[i]  =  u[i] + dt()/24.0*(  9.0* fCorrOld[i] 
                                              + 19.0* rhs.f[i](t,u) 
                                              -  5.0* rhs.f[i](t-dt(), um1 )
                                              +  1.0* rhs.f[i](t-(2*dt()), um2 )  );  



              fCorr[i] = rhs.f[i](t+dt() , uCorr);

              error[i] = fabs(uCorr[i]-uCorrOld[i]);

              fCorrOld[i] = fCorr[i] ;        

              iter++ ;    
           }
            
            up1[i] = uCorr[i] ; //- CORRECTOR  u(t+1) 
          } 
            
//---------- next step 
//
            for(auto i=0 ; i< u.size(); i++)
            {  
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
inline void AdamsMoulton4thSolver<Type>::solve() noexcept 
{
     std::cout << "Running Adams bashforth (4 step) - Adams Moulton (3 step) Solver" << std::endl;
      
         u.resize(u0().size());
         um1.resize(u0().size());
         um2.resize(u0().size());
         um3.resize(u0().size());
         up1.resize(u0().size());
         
         uPred.resize(u0().size());
         uCorr.resize(u0().size());
         uCorrOld.resize(u0().size());
         
         fPred.resize(u0().size());
         fCorr.resize(u0().size());
         fCorrOld.resize(u0().size());    // using to iterate 
         
         K1.resize(u0().size());
         K2.resize(u0().size());
         K3.resize(u0().size());
         K4.resize(u0().size());
         
         error.resize(u0().size());

         t = t0();

         std::cout << t << ' ' ;
         for(auto i=0 ; i < u.size() ; i++)
         {
            um3[i] = u0()[i] ;                  // set initial value 
            std::cout << um3[i] << ' ' ;                // and write them to file
         }
         std::cout << std::endl;
         
         
         
         for(auto i=0 ; i< u.size() ; i++)
         {  
            K1[i] = rhs.f[i]( t          ,  um3 );
            K2[i] = rhs.f[i]( t + dt()/2 ,  um3 + dt()/2 * K1[i] );
            K3[i] = rhs.f[i]( t + dt()/2 ,  um3 + dt()/2 * K2[i] );
            K4[i] = rhs.f[i]( t + dt()   ,  um3 + dt()   * K3[i] );
                                      // compute first point (start - up  Adams Bashforth )
            um2[i] =  um3[i] + dt()/6 *(K1[i] + 2*K2[i] + 2*K3[i] + K4[i]) ;  // Rk 4th order 
         }
            

         t+=dt();   
            
         std::cout << t << ' ' ;
         for(auto i=0 ; i < u.size() ; i++)
            std::cout << um2[i] << ' ' ;                // and write them to file
         std::cout << std::endl;
    
         for(auto i=0 ; i< u.size() ; i++)
         {  
            K1[i] = rhs.f[i]( t          ,  um2 );
            K2[i] = rhs.f[i]( t + dt()/2 ,  um2 + dt()/2*K1[i] );
            K3[i] = rhs.f[i]( t + dt()/2 ,  um2 + dt()/2*K2[i] );
            K4[i] = rhs.f[i]( t + dt()   ,  um2 + dt()  *K3[i] );
            
                                      // compute second point (start - up  Adams Bashforth )
            um1[i] =  um2[i] + dt()/6 *(K1[i] + 2*K2[i] +2*K3[i] + K4[i] ) ;  // Rk 4th order 
         }

         t+=dt();
         
         std::cout << t << ' ' ;
         for(auto i=0 ; i < u.size() ; i++)
            std::cout << um1[i] << ' ' ;                // and write them to file
         std::cout << std::endl;
         
         for(auto i=0; i < u.size() ; i++)
         {
            K1[i] = rhs.f[i]( t            , um1 );
            K2[i] = rhs.f[i]( t +   dt()/2 , um1+ dt()/2 *K1[i] );
            K3[i] = rhs.f[i]( t +   dt()/2 , um1+ dt()/2 *K2[i] );
            K4[i] = rhs.f[i]( t +   dt()   , um1+ dt()   *K3[i] );
                                                      // compute third point
            u[i]   = um1[i] + dt()/6 *(K1[i] + 2*K2[i] + 2*K3[i] + K4[i]);  // Rk 4th order 
         }

     /// @MAIN LOOP
     //
         for(t = t+dt() ; t < tf() ; t += dt() )
         {
            std::cout << t << "  " ;
            for(auto i=0 ; i < u.size() ; i++)
                  std::cout << u[i] << " " ;
            std::cout << std::endl; 


//------------ PREDICTOR      ADAMS-BASHFORTH 4TH ORDER (4 step)
//
            for(auto i=0 ; i< u.size(); i++)
            { 
                  uPred[i] = u[i] + dt()/24.0 *( 55.0 * rhs.f[i](t , u) 
                                               - 59.0 * rhs.f[i](t-dt() , um1) 
                                               + 37.0 * rhs.f[i](t-(2*dt()), um2 ) 
                                               -  9.0 * rhs.f[i](t-(3*dt()), um3 ) );
            

                  fPred[i] = rhs.f[i](t+dt() , uPred);     
            }

//------------ CORRECTOR     ADAMS-MOULTON 4TH ORDER (3 step)          
            
            for(auto i=0; i< u.size() ; i++){
                  error[i] = 1.0 ;        
                  fCorrOld[i] = fPred[i] ;
            }

            std::size_t iter = 0 ;
           
          for(auto i=0 ; i < u.size() ; i++)
          {
           while( error[i] >= toll )
           {
              uCorrOld[i] = u[i] + dt()/24.0*(   9.0* fCorrOld[i] 
                                              + 19.0* rhs.f[i](t,u) 
                                              -  5.0* rhs.f[i](t-dt(), um1 )
                                              +  1.0* rhs.f[i](t-(2*dt()), um2 )  );  
                  
              fCorrOld[i] = rhs.f[i](t+dt(), uCorrOld );

              uCorr[i]  =  u[i] + dt()/24.0*(  9.0* fCorrOld[i] 
                                              + 19.0* rhs.f[i](t,u) 
                                              -  5.0* rhs.f[i](t-dt(), um1 )
                                              +  1.0* rhs.f[i](t-(2*dt()), um2 )  );  



              fCorr[i] = rhs.f[i](t+dt() , uCorr);

              error[i] = fabs(uCorr[i]-uCorrOld[i]);

              fCorrOld[i] = fCorr[i] ;        

              iter++ ;    
           }
            
            up1[i] = uCorr[i] ; //- CORRECTOR  u(t+1) 
          } 
            
//---------- next step 
//
            for(auto i=0 ; i< u.size(); i++)
            {  
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
