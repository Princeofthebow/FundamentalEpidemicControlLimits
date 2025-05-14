The `Matlab` code in `SI_Model_ValueIteration.m` finds the exact optimal policies for the considered SI model for a specific set of parameters. In order to obtain the optimal policies the user can run the script using the desired set of parameters.
The code in `SI_Gillespie.m` simulates using the Stochastic simulation algorithm for considered Markov jump process version of the SI epidemic.


The `C++` code in the files `SI_Bursts_PolicyComputation` and `SI_PolicyComputation` finds the exact optimal policies for the considered SI model with and without the bursts over a range of parameters of interest.
The code has been written and tested in a Linux environment using the `g++` compiler for `C++`.
In order to compile the file the user will need `g++` as well as the `openMp` library. In order to execute the code over a specific set of parameters the user has to edit the file and modify the appropriate parameter grids and variables of interest.


In order to compile the code the following commands has to be executed in a terminal in the folder containing the file of interest.

    ```console
    foo@bar:~$ g++ -fopenmp -O3 -o filetorun.o hello.cpp 
    foo@bar:~$ ./filetorun.o
    ```

*N.B.* The computation time can be long in certain parameter configurations and for large parameter grids. The default parameters are set to obtain a solution in a relatively short amount of time on a modern computer. However, in order to obtain the optimal policies for large populations and a wide range of parameters, it is advised to:


Once the computation has finished each of the file will output the results to several txt files which can processed with the `Matlab` scripts available in the appropriate folder in `scripts_plots`.
