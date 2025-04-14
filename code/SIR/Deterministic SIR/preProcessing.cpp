// This code is published under the Eclipse Public License
// Authors: Pierre Martinon, Benjamin Heymann, Frederic Bonnans
// Inria Saclay and CMAP Ecole Polytechnique
// 2017

/**
  *\fn int preProcessing(const int dim_constants,
                         const vector<double>& constants,
                         vector<double>& starting_point,
                         int& starting_mode)
  * \param dimConstants  : number of constants
  * \param constants     : array of constants
  * \param starting_point: array of initial conditions for trajectory simulation
  * \param starting_mode: initial discrete mode for trajectory simulation
  *
  * Called once before starting the optimization.
  * Can be used to set initial conditions for the trajectory simulation.
  * Can also be used for instance to read external data in global (extern) variables
  */

#include "header_preProcessing"
{
    return 0;
}
