// This code is published under the Eclipse Public License
// Authors: Daphne Giorgi, Benjamin Heymann, Jinyan Liu, Pierre Martinon, Olivier Tissot
// Inria Saclay and Cmap Ecole Polytechnique
// 2014-2017


// Function for the state admissibility
// Input :
// time : current time (t)
// state : vector of state variables (x)
// mode : current mode of the system (i)
// constants : vector of constants
// dim_constant : dimension of the vector constants
// Output :
// true if the state is admissible
// false if it is not
#include "header_checkAdmissibleState"
{

  return true;
}


// Function for the (control,state) admissibility
// Input :
// time : current time (t)
// state : vector of state variables (x)
// control: vector of control variables (u)
// mode : current mode of the system (i)
// constants : vector of constants
// dim_constant : dimension of the vector constants
// Output :
// true if the (control,state) pair is admissible
// false if it is not
#include "header_checkAdmissibleControlState"
{
  double s = state[0];
  double i = state[1];
  const double N = constants[0];

  if(s+i<=N){
  return true;}
	else{
  return false;}

}


// Function for the final state admissibility
// Input :
// time : current time (t)
// final_state : vector of state variables (x)
// control: vector of control variables (u)
// mode : current mode of the system (i)
// constants : vector of constants
// dim_constant : dimension of the vector constants
// Output :
// true if the (control,state) pair is admissible
// false if it is not
#include "header_checkAdmissibleFinalState"
{
  return true;
}
