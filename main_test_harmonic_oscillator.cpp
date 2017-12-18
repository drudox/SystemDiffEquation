# include <iostream>
# include <iomanip>
# include <string>
# include "rhsODEproblem.H" 
# include <cmath>
# include "Euler/ForwardEulerSolver.H"
# include "Euler/BackwardEulerSolver.H"
# include "RungeKutta/ModifiedEulerSolver.H"
# include "RungeKutta/HeunSolver.H"
# include "RungeKutta/RungeKutta4th/RungeKutta4thSolver.H"
# include "Multistep/AdamsMethods/AdamsBashforth/AdamsBashforth2ndSolver.H"

using namespace std;
using namespace mg::numeric::odesystem ;
 

/*------------------------------------------------------------------
 *
 *      Test System of IVP harmonic oscillator 
 *
 *
 *      @Marco Ghiani Dec 2017 Glasgow    
 *
 -----------------------------------------------------------------*/


auto exacFun = [](double t, double y) { return exp(-5*pow((t-1),2) ); };

auto numFun1 =[](const double t , std::valarray<double> y){return  y[1] ; } ;
auto numFun2 =[](const double t , std::valarray<double> y){return  -y[0] -0.15*y[1] ; } ;


int main()
{
   
   std::vector<std::function<const double (const double, const std::valarray<double>)>>  function ;
   
   std::vector<std::function<const double (const double, const double)>>  real ;
   
   real.push_back(exacFun);

   function.push_back(numFun1);
   function.push_back(numFun2);
  
   const double t0 = 0.0;
   const double tf = 100.0 ;
   const double dt = 0.01;
   const std::valarray<double> u0 = {1.0,0.0} ;
      
   string fname = "analitical_1.out" ; 
    
    
   rhsODEProblem<double> p1(function, t0, tf , dt, u0 );
         
   ForwardEulerSolver<double> feuler1(p1) ;
   feuler1.solve("fwdEulerSystem_ocillator.out") ;
   
   BackwardEulerSolver<double> beuler1(p1) ;
   beuler1.solve("bwdEulerSystem_oscillator.out") ;
   




   ModifiedEulerSolver<double> meuler1(p1) ;
   meuler1.solve("ModEulerSystem_oscillator.out") ;

   HeunSolver<double> heun1(p1);   
   heun1.solve("HeunSystem_oscillator.out"); 

   RungeKutta4Solver<double> rk4_1(p1);
   rk4_1.solve("RK4System_oscillator.out");   

   
   
  return 0;    

} 


