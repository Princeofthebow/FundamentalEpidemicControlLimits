// This code is published under the Eclipse Public License
// Authors: Daphne Giorgi, Benjamin Heymann, Jinyan Liu, Pierre Martinon, Olivier Tissot
// Inria Saclay and Cmap Ecole Polytechnique
// 2014-2017

// General dynamics
// dy/dt = drift(t,y,u)dt + volatility(y,u)dWt where Wt is the standard Brownian motion

// Function for the drift (deterministic dynamics)
// Input :
// time : current time t
// initial_time : t0
// final_time : tf
// state : vector of state variables x
// control : vector of control variables u
// mode : mode of the system i
// constants : vector of constants
// dim_constant : dimension of the vector constants
// Output :
// state_dynamics : drift f(t,x,u,i) ie deterministic dynamics
#include "header_drift"
{

  //Example: double integrator, state = (x1 x2), control = (u) 
  //\dot x1 = x2
  //\dot x2 = u

  double s = state[0];
  double i = state[1];
  double u = control[0];
  const double N = constants[0];
  const double gamma = constants[5];
  const double beta = constants[6];

  state_dynamics[0] = -beta*s*i*(1-u);
  state_dynamics[1] = +beta*s*i*(1-u)-gamma*i;


}


// Function for the volatility (stochastic dynamics)
// Input :
// time : current time t
// initial_time : t0
// final_time : tf
// state : vector of state variables x
// control : vector of control variables u
// mode : mode of the system i
// constants : vector of constants
// dim_constant : dimension of the vector constants
// Output :
// volatility_dynamics : vector giving the volatility expression of the volatility
// Remember that this is a matrix of dimension dim_state x dim_brownian and you have to fill every coefficient.
#include "header_volatility"
{
}
