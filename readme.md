# Documentation on the code related to the paper *Fundamental limits for taming infectious disease*


## Table of Contents
- [Project description](#introduction)
- [Directory structure](#directory-structure)

## Project description

**Contributors**

This research was conceived and led by Giovanni Pugliese Carratelli and Ioannis Lestas with the contribution of Kris Varun Parag and Xiaodong Cheng.

**Summary**

This repository provides the code, some examples and the plots for the paper *Fundamental limits for taming infectious disease* which investigates how during an epidemic prevalence signals can be ineffective to tame the spread of a disease.

We consider the problem of finding transmission mitigation measures $u^\ast$ that minimise the social/economic costs $g_c(i,u)$ due to $i$ infected and intervention $u$ over an indefinite amount of time.

![Fundamental limits on taming epidemics](Diagram.png)

The paper demonstrates that during an epidemic increasing infections are not associated with increasingly restrictive measures for large class of relevant costs and policies are predominately constant. To this end we leverage stochastic jump models that accurately describe the evolution of the epidemics in that only rates can be controlled rather than individual events. We then make use of the Bellman Equation and we develop computational and analytical tools and we show that for a large class of models and cost functions prevalence signals are not effective.


This repository contains the scripts and data files used to generate the results presented in the paper *Fundamental limits for taming infectious disease* and its supplementary materials. Below is an overview of the directory structure and a list of computer source, and data processing scripts.

## Directory Structure and summary of implementatino details

There are two main folders, `code/` and `plots/`. `code` contains all the scripts used to generate the data files. These include ```C++``` code as well as ```Matlab``` code. `script_plots` contains several data processing files and will contain the code and material that is needed to reconstruct the figures of the paper and the supplementary material.
Specific requirement for the codebase are provided in the folders outlined below but in summary:

- The ```C++``` code make use of the `openmp` `g++` compiler directive and in order to obtain the the optimal policies for large populations and a large parameter ranges it is advised to run it in High Perfomance Computing environment. 
- The ```Matlab``` code does not require any specific toolbox.

---

The content of `code/` and as well as a short description of the content is provided below.
  - **SIR/**: Contains the ```Matlab``` code and functions related to computation of the optimal policies as well as the stochastic simulations for the time-evolution of the Markov Jump process version of the *Stochastic Infectious Recovered* model we consider.
  - **SEIR/**: Contains the ```Matlab``` code and functions related to computation of the optimal policies as well as the stochastic simulations for the time-evolution of the Markov Jump process version of the *Stochastic Exposed Infectious Recovered* models we consider.
  - **SI/**: Contains the ```C++``` and  ```Matlab``` code and functions related to the computation of the optimal policies as well as the stochastic simulations for the time-evolution of the Markov Jump Process version of the *Stochastic Infectious* models we consider.

<details>
<summary> Content of `code/SI` </summary>

- **`readme.md`**: Readme file providing the requirements to run the code and guidance on how to run the files in the folder and/or change parameter over which to do computations. The file also provides guidance how to process the resulting data.

- **`SI_Model_ValueIteration.m`**: `Matlab` function that will compute the optimal policy for the SI model for a specific parameter configuration single 

- **`SI_PolicyComputation.cpp`**: `C++` file that will compute the optimal policy for the SI model for a specific range of parameters of interest 

- **`SI_Bursts_PolicyComputation.cpp`**: `C++` file that will compute the optimal policy for the SI model when considering infection bursts (i.e. when an infection event can lead to a random number of infected individuals) for a specific range of parameters of interest

- **`SI_Gillespie.m`**: `Matlab` script implementing the Stochastic Simulation Algorithm in order to obtain the trajectories of an SI epidemic. The function runs multiple times in parallel for specific set of parameters and can be adapted to incorporate parameter variations.
</details>

<details>
<summary> Content of `code/SIR` </summary>

- **`readme.md`**: Readme file providing the requirements to run the code and guidance on how to run the files in the folder and/or change parameter over which to do computations. The file also provides guidance how to process the resulting data.

- **`SIR_PolicyComputation.m`**: `Matlab` script that will compute the optimal policy for the SIR model for a specific parameter configuration. 

- **`SIR_BurstPolicyComputation.m`**: `Matlab` script that will compute the optimal policy for the SIR model for a specific parameter configuration i the presence of bursts, i.e. when an infection event can lead to a random number of infected individuals. 

- **`SIR_Gillespie.m`**: `Matlab` script implementing the Stochastic Simulation Algorithm in order to obtain the trajectories of an SIR epidemic. The function runs multiple times in parallel for specific set of parameters and can be adapted to incorporate parameter variations.

</details>

<details>
<summary> Content of `code/SEIR` </summary>

- **`readme.md`**: Readme file providing the requirements to run the code and guidance on how to run the files in the folder and/or change parameter over which to do computations. The file also provides guidance how to process the resulting data.

- **`SEIR_PolicyComputation.m`**: `Matlab` function that will compute the optimal policy for the SEIR model for a specific parameter configuration. 

</details>
---

The folder `script_plots/` contains data processing files  and after the review process is finalised will contain sub-folders that are named after figures in the paper. The sub-folders will include the relevant figure, the data and the code that was used to generate the figure.