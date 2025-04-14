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

The paper demonstrates that during an epidemic increasing infections are not associated with increasingly restrictive measures for large class of relevant costs. To this end we leverage stochastic jump models that accurately describe the evolution of the epidemics in that only rates can be controlled rather than individual events. We then make use of the Bellman Equation and we develop computational and analytical tools and we show that for a large class of models and cost functions prevalence signals are not effective.


This repository contains the scripts and data files used to generate the results presented in the paper *Fundamental limits for taming infectious disease* and its supplementary materials. Below is an overview of the directory structure and a list of computer source, and data processing scripts.

## Directory Structure

There are two main folders, `code/` and `plots/`. `code` contains all the scripts used to generate the data files. These include ```C++``` code as well as ```Matlab``` code. `script_plots` contains several data processing files and will contain the code and material that is needed to reconstruct the figures of the paper and the supplementary material.

---

The content of `code/` and as well as a short description of the content is provided below.
  - **SIR/**: Contains the ```Matlab``` code and functions related to computation of the optimal policies as well as the stochastic simulations for the time-evolution of the Markov Jump process version of the *Stochastic Infectious Recovered* model we consider.
  - **SEIR/**: Contains the Matlab code and functions related to computation of the optimal policies as well as the stochastic simulations for the time-evolution of the Markov Jump process version of the *Stochastic Exposed Infectious Recovered* models we consider.
  - **SI/**: Contains the ```C++``` code and functions related to the computation of the optimal policies as well as the stochastic simulations for the time-evolution of the Markov Jump Process version of the *Stochastic Infectious* models we consider.

<details>
<summary> Content of `code/SI` </summary>

- **`readme.md`**: Readme file providing guidance on how to run the files in the folder and/or change parameter over which to do computations. The file also provides guidance how to process the resulting data.

- **`SI_Model_ValueIteration.m`**: `Matlab` function that will compute the optimal policy for the SI model for a specific parameter configuration single 

- **`SI_PolicyComputation.cpp`**: `C++` file that will compute the optimal policy for the SI model for a specific range of parameters of interest 

- **`SI_Bursts_PolicyComputation.cpp`**: `C++` file that will compute the optimal policy for the SI model when considering infection bursts for a specific range of parameters of interest 
</details>

<details>
<summary> Content of `code/SIR` </summary>

- **`readme.md`**: Readme file providing guidance on how to run the files in the folder and/or change parameter over which to do computations. The file also provides guidance how to process the resulting data.

- **`SIR_Model_ValueIteration.m`**: `Matlab` function that will compute the optimal policy for the SIR model for a specific parameter configuration. 

</details>

<details>
<summary> Content of `code/SEIR` </summary>

- **`readme.md`**: Readme file providing guidance on how to run the files in the folder and/or change parameter over which to do computations. The file also provides guidance how to process the resulting data.

- **`SEIR_PolicyComputation.m`**: `Matlab` function that will compute the optimal policy for the SEIR model for a specific parameter configuration. 

</details>
---

The folder `script_plots/` contains data processing files  and after the review process is finalised will contain sub-folders that are named after figures in the paper. The sub-folders will include the relevant figure, the data and the code that was used to generate the figure.