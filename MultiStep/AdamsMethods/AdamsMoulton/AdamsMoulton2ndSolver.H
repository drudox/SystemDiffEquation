# ifndef __ADAMS_MOULTON_2STEP_PC_SOLVER_H__
# define __ADAMS_MOULTON_2STEP_PC_SOLVER_H__


# include "../AdamsMethods.H"
# include "../../../rhsODEproblem.H"


namespace mg {
                namespace numeric {
                                     namespace odesystem {


/**------------------------------------------------------------------------------------
 *
 *\ @brief \class AdamsMoulton 2nd order Solver
 * 
 *  Perform linear multi-step (2) scheme using PREDICTOR (Adams-bashforth 2 step) 
 *                                           & CORRECTOR (Adams-Moulton (imp) 2 step )  
 * 
 *
 *\ @Marco Ghiani Dec 2017, Glasgow UK
 *
 *
 ------------------------------------------------------------------------------------*/



template<typename Type = double>
class AdamsMoulton2ndSolver :
                                 public AdamsMethods<Type> 
{

//---      
    public:  

      AdamsMoulton2ndSolver(const rhsODEProblem<Type> & that) noexcept :
                                                                          AdamsMethods<Type>{that} 
                  {}
      
      virtual ~AdamsMoulton2ndSolver() = default ;

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
      
      using AdamsMethods<Type>::error;


};

//------------------  Implementation (to be put into .cpp file)   -----------------  //


template<typename Type>
inline void AdamsMoulton2ndSolver<Type>::solve(const std::string filename) {
      
    std::ofstream f(filename, std::ios::out );
      
    if(!f)
    {     
       throw std::runtime_error(">> Error opening File in ADAMS BASHFORTH-MOULTON  solver <<" );
    }
    else
    {
         std::cout << "Running Adams Bashforth (2 step) - Adams Moulton (1 step) Solver" << std::endl;
      
         u.resize(u0().size());
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

         f << t << ' ' ;
         for(auto i=0 ; i < u.size() ; i++)
         {
            um1[i] = u0()[i] ;                  // set initial value 
            f << um1[i] << ' ' ;                // and write there to file
         }
         f << std::endl;
         
            
         for(auto i=0 ; i< u.size() ; i++)
         {  
            K1[i] = rhs.f[i]( t        ,  um1 );
            K2[i] = rhs.f[i]( t + dt() ,  um1 + dt()*K1[i] );
                                      // compute first point (start - up  Adams Bashforth )
            u[i] =  um1[i] + dt()/2 *(K1[i] + K2[i]) ;  // Rk 2nd order 
         }


         for(t+=dt() ; t <= tf() ; t += dt() )
         {
            f << t << "  " ;
            for(auto i=0 ; i < u.size() ; i++)
                  f << u[i] << " " ;
            f << std::endl; 


//------------ PREDICTOR  ADAMS BASHFORTH 2th ORDER (2 STEP)
//
            for(auto i=0 ; i< u.size(); i++)
            { 
                  uPred[i] = u[i] + dt()/2 *( 3.0 * rhs.f[i](t , u) 
                                             -1.0 * rhs.f[i](t-dt() , um1) );
            

                  fPred[i] = rhs.f[i](t+dt() , uPred);     
            }

//------------ CORRECTOR  ADAMS MOULTON 2th ORDER (1 STEP)   became the TRAPEZOIDAL (C-N scheme)        
            
            for(auto i=0; i< u.size() ; i++){
                  error[i] = 1.0 ;        
                  fCorrOld[i] = fPred[i] ;
            }

            std::size_t iter = 0 ;
           
          for(auto i=0 ; i < u.size() ; i++)
          {
           while( error[i] >= toll )
           {
                  uCorrOld[i] = u[i] + dt()/2.0*( 1.0 * fCorrOld[i] 
                                                + 1.0 * rhs.f[i](t,u) ) ;  
                  
                  fCorrOld[i] = rhs.f[i](t+dt(), uCorrOld );

                  uCorr[i]  =  u[i] + dt()/2.0*( 1.0 * fCorrOld[i] 
                                                +1.0 * rhs.f[i](t,u) ); 
               
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
               um1[i] =   u[i];
                 u[i] = up1[i];
            } 
         }  
         std::cout << "... Done " << std::endl;  
         f.close();
    } 
}


template<typename Type>
inline void AdamsMoulton2ndSolver<Type>::solve() noexcept 
{
 
     std::cout << "Running Adams bashforth-Moulton (2 step) Solver" << std::endl;
      
         u.resize(u0().size());
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
            um1[i] = u0()[i] ;                  // set initial value 
            std::cout << um1[i] << ' ' ;                // and write there to file
         }
         std::cout << std::endl;
         
            
         for(auto i=0 ; i< u.size() ; i++)
         {  
            K1[i] = rhs.f[i]( t        ,  um1 );
            K2[i] = rhs.f[i]( t + dt() ,  um1 + dt()*K1[i] );
                                      // compute first point (start - up  Adams Bashforth )
            u[i] =  um1[i] + dt()/2 *(K1[i] + K2[i]) ;  // Rk 2nd order 
         }

      /// @Main LOOP
         for(t+=dt() ; t <= tf() ; t += dt() )
         {
            std::cout << t << "  " ;
            for(auto i=0 ; i < u.size() ; i++)
                  std::cout << u[i] << " " ;
            std::cout << std::endl; 


//------------ PREDICTOR  ADAMS BASHFORTH 2th ORDER (2 STEP)
//
            for(auto i=0 ; i< u.size(); i++)
            { 
                  uPred[i] = u[i] + dt()/2 *( 3.0 * rhs.f[i](t , u) 
                                             -1.0 * rhs.f[i](t-dt() , um1) );
            

                  fPred[i] = rhs.f[i](t+dt() , uPred);     
            }

//------------ CORRECTOR  ADAMS MOULTON 2th ORDER (1 STEP)   became the TRAPEZOIDAL (C-N scheme)        
            
            for(auto i=0; i< u.size() ; i++){
                  error[i] = 1.0 ;        
                  fCorrOld[i] = fPred[i] ;
            }

            std::size_t iter = 0 ;
           
          for(auto i=0 ; i < u.size() ; i++)
          {
           while( error[i] >= toll )
           {
                  uCorrOld[i] = u[i] + dt()/2.0*( 1.0 * fCorrOld[i] 
                                                + 1.0 * rhs.f[i](t,u) ) ;  
                  
                  fCorrOld[i] = rhs.f[i](t+dt(), uCorrOld );

                  uCorr[i]  =  u[i] + dt()/2.0*( 1.0 * fCorrOld[i] 
                                                +1.0 * rhs.f[i](t,u) ); 
               
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
