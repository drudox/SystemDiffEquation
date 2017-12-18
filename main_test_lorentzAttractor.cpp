# include <iostream>
# include <iomanip>
# include <string>
# include "rhsODEproblem.H" 
# include <cmath>
# include "ForwardEulerSolver.H"
# include "BackwardEulerSolver.H"
# include "ModifiedEulerSolver.H"
# include "HeunSolver.H"
# include "RungeKutta4Solver.H"
# include "LeapFrogSolver.H"
# include "AdamsMoulton2ndSolver.H"
# include "AdamsBashforth2ndSolver.H"


using namespace std;
using namespace mg::numeric::odesystem ;
 

/*------------------------------------------------------------------
 *
 *      Test System of IVP lorentz attractor 
 *
 *
 *      @Marco Ghiani Dec 2017 Glasgow    
 *
 -----------------------------------------------------------------*/


auto exacFun = [](double t, double y) { return exp(-5*pow((t-1),2) ); };

auto numFun1 =[](const double t , std::valarray<double> y){return  10.0 * (y[1] - y[0]) ; } ;
auto numFun2 =[](const double t , std::valarray<double> y){return  28.0 * y[0] - y[1] - y[0] * y[2] ; } ;
auto numFun3 =[](const double t , std::valarray<double> y){return -8.0/3.0 * y[2] + y[0]* y[1] ; } ;


int main(){
   
   std::vector<std::function<const double (const double, const std::valarray<double>)>>  function ;
   
   std::vector<std::function<const double (const double, const double)>>  real ;
   
   real.push_back(exacFun);

   function.push_back(numFun1);
   function.push_back(numFun2);
   function.push_back(numFun3);
  
   const double t0 = 0.0;
   const double tf = 100.0 ;
   const double dt = 0.01;
   const std::valarray<double> u0 = {1.0,0.0,0.0} ;
      
   string fname = "analitical_1.out" ; 
    
    
   rhsODEProblem<double> p1(function, t0, tf , dt, u0 );
         
   ForwardEulerSolver<double> feuler1(p1) ;
   feuler1.solve("fwdEulerSystem_lorentz.out") ;
   /*
   BackwardEulerSolver<double> beuler1(p1) ;
   beuler1.solve("bwdEulerSystem_pp.out") ;
   */

   AdamsBashforth2ndSolver<double> ab2(p1);
   ab2.solve("AdamsBashforth2ndSystem_lorentz.out");

   AdamsMoulton2ndSolver<double> am2(p1);
   am2.solve("AdamsMoulton2ndSystem_lorentz.out");


   ModifiedEulerSolver<double> meuler1(p1) ;
   meuler1.solve("ModEulerSystem_lorentz.out") ;

   HeunSolver<double> heun1(p1);   
   heun1.solve("HeunSystem_lorentz.out"); 

   RungeKutta4Solver<double> rk4_1(p1);
   rk4_1.solve("RK4System_lorentz.out");   

   LeapFrogSolver<double> leapFrog(p1);
   leapFrog.solve("LeapFrog_lorentz.out");
   
  return 0;    
}
