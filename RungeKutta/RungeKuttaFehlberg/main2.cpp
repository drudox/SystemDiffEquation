# include <iostream>
# include <iomanip>
# include <string>
# include "../../rhsODEproblem.H" 
# include <cmath>

#include "../RungeKutta.H"
# include "RungeKuttaFehlberg45Solver.H"
# include "../RungeKutta4th/RungeKutta4thSolver.H"



using namespace std;
using namespace mg::numeric::odesystem ;
 



//auto numFun1 =[](const double t , std::vector<double> y){return  1.0 * y[0] - 1.0* y[0] * y[1] ; } ;
//auto numFun2 =[](const double t , std::vector<double> y){return -0.5 * y[1] + 0.5* y[0] * y[1] ; } ;
auto numFun = [](double t, std::valarray<double> u) { return -10*(t-1)*u[0]; } ;

auto exacFun = [](double t, double y) { return exp(-5*pow((t-1),2) ); };


int main(){
   
   std::vector<std::function<const double (const double, const std::valarray<double>)>>  function ;
   
   std::vector<std::function<const double (const double, const double)>>  real ;
   
   real.push_back(exacFun);

   function.push_back(numFun);
 //  function.push_back(numFun2);
  
   const double t0 = 0.0;
   const double tf = 2 ;
   const double dt = 0.04;
 //  const std::vector<double> u0 = {1 ,0.1} ;
   const std::valarray<double> u0 = {exp(-5.0)} ;  
   string fname = "analitical_1.out" ; 
    
    
   rhsODEProblem<double> p1(function, t0, tf , dt, u0 );
         
   
 

   RungeKutta4Solver<double> rk4(p1); 
   rk4.solve("rk4_1.out");

   RungeKuttaFehlberg45Solver<double> rkf1(p1); 
   rkf1.solve("rkf_1.out");
      
   
  return 0;    
}
