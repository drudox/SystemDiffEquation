# ifndef __ADAMS_MOULTON_3STEP_PC_SOLVER_H__
# define __ADAMS_MOULTON_3STEP_PC_SOLVER_H__


# include "../AdamsMethods.H"
# include "../../../rhsODEproblem.H"


namespace mg {
                namespace numeric {
                                     namespace odesystem {


/**------------------------------------------------------------------------------------
 *
 *\ @brief \class AdamsMoulton3thSolver
 * 
 *  Perform linear multi-step scheme using PREDICTOR (Adams-bashforth 3 step : 3th ord)  
 *                                       & CORRECTOR (Adams-Moulton (imp) 2 step : 3th ord )  
 * 
 *
 *\ @Marco Ghiani Dec 2017, Glasgow UK
 *
 *
 ------------------------------------------------------------------------------------*/



template<typename Type = double>
class AdamsMoulton3thSolver :
                                 public AdamsMethods<Type> 
{

//---      
    public:  

      AdamsMoulton3thSolver(const rhsODEProblem<Type> & that) noexcept :
                                                                          AdamsMethods<Type>{that} 
                  {}
      
      virtual ~AdamsMoulton3thSolver() = default ;

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
inline void AdamsMoulton3thSolver<Type>::solve(const std::string filename) {
      
    std::ofstream f(filename, std::ios::out );
      
    if(!f)
    {     
       throw std::runtime_error(">> Error opening File in ADAMS BASHFORTH-MOULTON  solver <<" );
    }
    else
    {
         std::cout << "Running Adams Bashforth (3 step) - Adams Moulton (2 step) Solver" << std::endl;
      
         u.resize(u0().size());
         um2.resize(u0().size());
         um1.resize(u0().size());
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
            um2[i] = u0()[i] ;                  // set initial value 
            f << um2[i] << ' ' ;                // and write there to file
         }
         f << std::endl;
         
            
         for(auto i=0 ; i< u.size() ; i++)
         {  
            K1[i] = rhs.f[i]( t          ,  um2 );
            K2[i] = rhs.f[i]( t + dt()/2 ,  um2 + dt()/2*K1[i] );
            K3[i] = rhs.f[i]( t + dt()/2 ,  um2 + dt()/2*K2[i] );
            K4[i] = rhs.f[i]( t + dt()   ,  um2 + dt()  *K3[i] );
 
            // compute first point (start - up  Adams Bashforth )
            um1[i] =  um2[i] + dt()/6 *(K1[i] + 2*K2[i] + 2*K3[i] + K4[i] ) ;  // Rk 4th order 
         }


         t += dt();

         f << t << ' ' ;
         for(auto i=0 ; i < u.size() ; i++)
            f << um1[i] << ' ' ;                // and write there to file
         f << std::endl;
         
            
         for(auto i=0 ; i< u.size() ; i++)
         {  
            K1[i] = rhs.f[i]( t          ,  um1 );
            K2[i] = rhs.f[i]( t + dt()/2 ,  um1 + dt()/2 *K1[i] );
            K3[i] = rhs.f[i]( t + dt()/2 ,  um1 + dt()/2 *K2[i] );
            K4[i] = rhs.f[i]( t + dt()   ,  um1 + dt()   *K3[i] );
        
        
        // compute second point (start - up  Adams Bashforth )
            u[i] =  um1[i] + dt()/6 *(K1[i] + 2*K2[i] + 2*K3[i] + K4[i]) ;  // Rk 4th order 
         }

      /// @Main LOOP
      //
         for(t= t+dt() ; t <= tf() ; t += dt() )
         {
            f << t << "  " ;
            for(auto i=0 ; i < u.size() ; i++)
                  f << u[i] << " " ;
            f << std::endl; 


//------------ PREDICTOR  ADAMS BASHFORTH 3th ORDER (3 STEP)
//
            for(auto i=0 ; i< u.size(); i++)
            { 
                  uPred[i] = u[i] + dt()/12 *(  23.0 * rhs.f[i](t , u) 
                                              - 16.0 * rhs.f[i](t-dt()   , um1) 
                                              +  5.0 * rhs.f[i](t-(2*dt()), um2) )  ;
            

                  fPred[i] = rhs.f[i](t+dt() , uPred);     
            }

//------------ CORRECTOR  ADAMS MOULTON 3th ORDER (2 STEP)          
            
            for(auto i=0; i< u.size() ; i++){
                  error[i] = 1.0 ;        
                  fCorrOld[i] = fPred[i] ;
            }

            std::size_t iter = 0 ;
           
          for(auto i=0 ; i < u.size() ; i++)
          {
           while( error[i] >= toll )
           {
                  uCorrOld[i] = u[i] + dt()/12*(5* fCorrOld[i] + 8 * rhs.f[i](t,u) - rhs.f[i](t-dt(), um1 ) );  
                  
                  fCorrOld[i] = rhs.f[i](t+dt(), uCorrOld );

                  uCorr[i]  =  u[i] + dt()/12*(5* fCorrOld[i] + 8 * rhs.f[i](t,u) - rhs.f[i](t-dt(), um1 ) ); 
               
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
inline void AdamsMoulton3thSolver<Type>::solve() noexcept 
{
    std::cout << "Running Adams bashforth-Moulton 3th order Solver" << std::endl;
      
         u.resize(u0().size());
         um2.resize(u0().size());
         um1.resize(u0().size());
         up1.resize(u0().size());
         
         uPred.resize(u0().size());
         uCorr.resize(u0().size());
         uCorrOld.resize(u0().size());
         
         fPred.resize(u0().size());
         fCorr.resize(u0().size());
         fCorrOld.resize(u0().size());    // using to iterate 
         
         K1.resize(u0().size());
         K2.resize(u0().size());
         
         error.resize(u0().size());
    
         t = t0();

         std::cout << t << ' ' ;
         for(auto i=0 ; i < u.size() ; i++)
         {
            um2[i] = u0()[i] ;                  // set initial value 
          std::cout << um2[i] << ' ' ;                // and write there to file
         }
         std::cout << std::endl;
         
            
         for(auto i=0 ; i< u.size() ; i++)
         {  
            K1[i] = rhs.f[i]( t        ,  um2 );
            K2[i] = rhs.f[i]( t + dt() ,  um2 + dt()*K1[i] );
                                      // compute first point (start - up  Adams Bashforth )
            um1[i] =  um2[i] + dt()/2 *(K1[i] + K2[i]) ;  // Rk 2nd order 
         }


         t += dt();

         std::cout << t << ' ' ;
         for(auto i=0 ; i < u.size() ; i++)
            std::cout << um1[i] << ' ' ;                // and write there to file
         std::cout << std::endl;
         
            
         for(auto i=0 ; i< u.size() ; i++)
         {  
            K1[i] = rhs.f[i]( t        ,  um1 );
            K2[i] = rhs.f[i]( t + dt() ,  um1 + dt()*K1[i] );
                                      // compute first point (start - up  Adams Bashforth )
            u[i] =  um1[i] + dt()/2 *(K1[i] + K2[i]) ;  // Rk 2nd order 
         }

      /// @Main LOOP
      //
         for(t= t+dt() ; t <= tf() ; t += dt() )
         {
            std::cout << t << "  " ;
            for(auto i=0 ; i < u.size() ; i++)
                  std::cout << u[i] << " " ;
            std::cout << std::endl; 


//------------ PREDICTOR  ADAMS BASHFORTH 3th ORDER (3 STEP)
//
            for(auto i=0 ; i< u.size(); i++)
            { 
                  uPred[i] = u[i] + dt()/12 *(  23.0 * rhs.f[i](t , u) 
                                              - 16.0 * rhs.f[i](t-dt()   , um1) 
                                              +  5.0 * rhs.f[i](t-(2*dt()), um2) )  ;
            

                  fPred[i] = rhs.f[i](t+dt() , uPred);     
            }

//------------ CORRECTOR  ADAMS MOULTON 3th ORDER (2 STEP)          
            
            for(auto i=0; i< u.size() ; i++){
                  error[i] = 1.0 ;        
                  fCorrOld[i] = fPred[i] ;
            }

            std::size_t iter = 0 ;
           
          for(auto i=0 ; i < u.size() ; i++)
          {
           while( error[i] >= toll )
           {
                  uCorrOld[i] = u[i] + dt()/12*(5* fCorrOld[i] + 8 * rhs.f[i](t,u) - rhs.f[i](t-dt(), um1 ) );  
                  
                  fCorrOld[i] = rhs.f[i](t+dt(), uCorrOld );

                  uCorr[i]  =  u[i] + dt()/12*(5* fCorrOld[i] + 8 * rhs.f[i](t,u) - rhs.f[i](t-dt(), um1 ) ); 
               
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
