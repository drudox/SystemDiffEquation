# ifndef __RHS_ODE_PROBLEM_H__
# define __RHS_ODE_PROBLEM_H__

# include <string>
# include <functional>
# include <fstream>
# include <vector>
# include <valarray>
# include <iostream>
//# include "Jacobian.H"


namespace mg {
               namespace numeric {
                                    namespace odesystem {

/*-----------------------------------------------------------------------
 *   @brief Class RightHandSide - ODE/SDE problem,
 *    
 *    --> Type : single precision (float) , double precision (double)
 *
 *    dy/dt = RHS 
 *
 *    @ Marco Ghiani  Dec 2017 Glasgow UK
 ------------------------------------------------------------------------*/


template <typename Type = double>
class rhsODEProblem {
  
//
//  
//--
   using analysisFunction = std::function<const Type(const Type, const std::valarray<Type>)>;
  
   public:
     class dfdx ;

     friend class dfdx ;
  
   public:
   
   

 //---  
   rhsODEProblem(const std::vector<std::function<const Type(const Type,const std::valarray<Type>)>> numfun ,
                  const std::vector<std::function<const Type(const Type,const std::valarray<Type>)>> exactfun ,
                  const Type,const Type,const Type,const std::valarray<Type>,const std::string ) noexcept ;
        
 //-- no analitical solution  
   rhsODEProblem(const std::vector<std::function<const Type(const Type,const std::valarray<Type>)>> numfun,
                 const Type, const Type, const Type, const std::valarray<Type> ) noexcept ;
 
 
      virtual ~rhsODEProblem() = default ;
      
      rhsODEProblem(const rhsODEProblem &) = default ;
      rhsODEProblem(rhsODEProblem&& ) = default ;
      rhsODEProblem& operator=(const rhsODEProblem&) = default ;
      rhsODEProblem& operator=(rhsODEProblem&& ) = default ;


      const std::vector<std::function<const Type(const Type,const std::valarray<Type>)>> numericalFunction ;
      
      const std::vector<std::function<const Type(const Type,const std::valarray<Type>)>> analiticalFunction ;

      
      const std::vector<analysisFunction>&  f = numericalFunction ;   

      
      const auto dfdt(std::size_t indx , const Type t , const std::valarray<Type> u) {
            return (f[indx](t, u+eps )-f[indx](t,u))/eps ;
      }   
      
     
      auto setRhs (const std::vector<
                   std::function<const Type(const Type,const std::valarray<Type>)>> numfun) noexcept
      {     
            
        for(auto i=0 ; i < numfun.size() ; i++)
        {
          numericalFunction.push_back(numfun.at(i)) ; 
        }  
      } 
      
      auto setExact(const std::vector<std::function<const Type(const Type,const std::valarray<Type>)>> exactfun) noexcept
      {
          for(auto i=0 ; i < exactfun.size(); i++)
          {
            analiticalFunction.push_back(exactfun.at(i));
          }  
      }
    
      auto solveExact() noexcept ;

      const Type t0() const noexcept { return _t0 ;} 
      const Type tf() const noexcept { return _tf ;}
      const Type dt() const noexcept { return _dt ;}
      const std::valarray<Type> u0() const noexcept { return _u0 ;}
      const std::string fname ()     const noexcept { return filename ;}
      
      

//---
   private:
     
     Type  _t0 ;  //! start time
     Type  _tf ;  //! final time
     Type  _dt ;  //! time-step
     std::valarray<Type> _u0 ;  //! initial value array (SysDE) 

     std::string filename ; 

     const Type eps  = 1e-12  ; 
      

};

template <typename Type>
rhsODEProblem<Type>::rhsODEProblem ( const std::vector< 
                                          std::function<const Type(const Type,const std::valarray<Type>)>> numfun ,
                                     const std::vector<
                                          std::function<const Type(const Type,const std::valarray<Type>)>> exactfun ,
                                     const Type Ti,const Type Tf,const Type Dt, 
                                     const std::valarray<Type> U0,
                                     const std::string fname 
                                   ) 
                                    noexcept : numericalFunction{numfun} ,
                                            analiticalFunction{exactfun} ,
                                                                 _t0{Ti} ,
                                                                 _tf{Tf} ,
                                                                 _dt{Dt} ,
                                                                 _u0{U0} ,    // vector copy constructor
                                                         filename{fname} 
{
      solveExact();
}


template<typename Type>
rhsODEProblem<Type>::rhsODEProblem ( const std::vector< 
                                      std::function<const Type(const Type,const std::valarray<Type>)>> numfun ,
                                     const Type Ti,const Type Tf,const Type Dt,const std::valarray<Type> U0           
                                   ) 
                                        noexcept : numericalFunction{numfun} , 
                                                                     _t0{Ti} ,
                                                                     _tf{Tf} ,
                                                                     _dt{Dt} ,
                                                                     _u0{U0}  
                   {}                                   

//- if exist ( and gives ) compute the 
//     numerical-exact solution 
//
template<typename Type>
auto rhsODEProblem<Type>::solveExact() noexcept {
   
   const Type Ns = ( _tf -_t0 )/ _dt ;
    

   std::ofstream fn( filename , std::ios::out);
   
   if(!fn)
   {
      std::cerr << "Error opening file " << filename << " in rhsODEProblem::solveExact()" << std::endl;      
      exit(-1);
   }
   else
   {
       //std::cout << Ns << std::endl ;   
      Type time = _t0 ;
      std::valarray<Type> yt;
      yt.resize(_u0.size()); ;
      
      fn << time << ' ' ; 
      for(auto i=0 ; i < yt.size() ; i++)
      {
              yt[i] = _u0[i] ;    
              fn << yt[i] << ' ' ;
      }        
      fn << std::endl ;   
       
       for(std::size_t i=0 ; i < Ns ; i++)
       { 
           time += _dt ; 
           fn << time << ' ' ;
           for(auto j=0 ; j < analiticalFunction.size() ; j++) {  
               yt[j] = analiticalFunction.at(j)(time , yt );
               fn  << yt[j] << ' ' ;
           } 
           fn << std::endl ;   
       }   
       fn.close();
    }
      
}

  }//ode 
 }//numeric
}//mg 
# endif
