The `Matlab` code in `SI_Model_ValueIteration.m` finds the exact optimal policies for the considered SI model for a specific set of parameters. in order to obtain the optimal policies the user can run the script using the desired set of parameters.

The `C++` code in the files `SI_Bursts_PolicyComputation` and `SI_PolicyComputation` finds the exact optimal policies for the considered SI model with and without the presence of bursts over a range of parameters of interest.
The code has been written and tested in a Linux environment using the `g++` compiler for `C++`.
In order to compile the file the user will need `g++` as well as the `openMp` library. In order to execute the code over a specific set of parameters the user has to edit the file and modify the appropriate parameter grids and variables of interest.


In order to compile the code the following commands has to be executed in the folder containing the file of interest.

``` g++ <filename.cpp> -fopenmp -O3```

once the compiling processes finishes the code can be run executing in the terminal

```./a.out ```

n.b. by default this will be the name of the executable object once the compilation has finished.

Once the computation has finished each of the file will output the results to several txt files which can processes with the `Matlab` scripts available in the appropriate folder in `scripts_plots`.
