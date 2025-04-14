#include <iostream>
#include <algorithm>
#include <cfloat>
#include <math.h>
#include <fstream>
#include <vector>
using namespace std;
using std::vector;

vector<double> linspace(double a, double b, int n)
{
    // This function creates a linearly space vectors
    vector<double> array;
    if ((n == 0) || (n == 1) || (a == b))
        array.push_back(b);
    else if (n > 1)
    {
        double step = (b - a) / (n - 1);
        int count = 0;
        while (count < n)
        {
            array.push_back(a + count * step);
            ++count;
        }
    }
    return array;
}

int main()
{
    // Instantiation of output streams to write to file the results
    fstream myfile;
    fstream myfileU;
    fstream myfileUidx;
    fstream myfileD;
    myfile.open("VI_B.txt", fstream::out);    // output stream for results on data
    myfileU.open("VI_U_B.txt", fstream::out); // output stream for saving the set of possible control actions
    myfileUidx.open("VI_Uidx_B.txt", fstream::out);
    myfileD.open("RunData.txt", fstream::out);

    // Control set variables
    const int noa = 3;     // # of control signals to select from
    const double umax = 1; // maximum control signal value
    vector<double> u = linspace(0, umax, noa);

    // System parameters
    const int N_i = 30;          // number of possible infected
    const int N = N_i + 1;       // number of nodes N in the markov chain --> N_i = N-1
    const double gamma = 0.32;   // recovery parameter
    const double dt = 1e-3;      // inverse of uniformisation factor
    const double invdt = 1 / dt; // uniformisation factor

    // Computation variables
    const double eps = 1e-6;     // Threshold for the value iteration algorithm
    const int c1gridsteps = 1;   // Number of steps of the grid of c1 costs
    const int c2gridsteps = 5;   // Number of steps of the c2 costs
    const int c3gridsteps = 1;   // Number of steps of the c4 costs
    const int c4gridsteps = 1;   // Number of steps of the c3 cost
    const int betagridsteps = 1; // Number of steps for the parameter beta
    const int pwrggridsteps = 1;
    const int qsupportsize = 1;
    const int qdists = 1;
    // double Pqs[qdists][qsupportsize] =  {{1,0,0},{0.75,0.2,0.05},{0.55,0.4,0.05},{0.55,0.3,0.15}};
    double Pqs[qdists][qsupportsize] = {{1}};

    vector<double> c1grid = {0};             // definition pf the grid for the costs c1
    vector<double> c2grid = linspace(0.1,10,c2gridsteps);       // definition of the grid for the costs c2
    vector<double> betagrid = {1.6 * gamma}; // definition of the grid for the parameter beta
    vector<double> c3grid = {1};             // definition of the grid for the costs c3
    vector<double> c4grid = {0};             // definition of the grid for the costs c4
    vector<double> pwrggrid = {1};

    const int pwrc2u = 1;    // power for the term z(u)
    const int pwrc3i = 1;    // power for the term c3i
    const double pwrc4u = 1; // power for the term c4u

    // Computation of total number of simulations obtained by taking the product of grid sizes
    const int N_sims = qdists * pwrggridsteps * c1gridsteps * c2gridsteps * betagridsteps * c3gridsteps * c4gridsteps;

    cout << "Will soon start to solve: " << N_sims << " Problems" << std::endl;
    // Instantiation of the indexes that are used in the various openmp threads
    int c1idx = 0;
    int c2idx = 0;
    int qidx = 0;
    int betaidx = 0;
    int c3idx = 0;
    int c4idx = 0;
    int pwrgidx = 0;

    // Instantiation of the cost and model variable that are used in the various openmp threads
    double c1;
    double c2;
    double beta;
    double c3;
    double c4;
    double pwrg;

    double ResultsJ[N_sims][N] = {0};
    double Resultsu[N_sims][N] = {0};
    int Resultsuidx[N_sims][N] = {0};

    // Instantiation of the result variables, these are both for the optimal total cost as well as the optimal policy
    double tempvar[noa] = {0}; // Temporary vector for VI minimisation
    double J[N] = {0};
    double u_opt[N] = {0};
    u_opt[0] = umax;
    double u_opt_idx[N] = {0};
    u_opt_idx[0] = noa; // in matlab indexes
    double delta = 0;
    double nonminterm = 0;
    double vold = 0; // Temp var for VI
    int min_idx = 0;
    double jmpterm = 0;

// The following directive runs the Value Iteration algorithm in parallel for each value
// of the considered cost/parameter in the grids dfined above
#pragma omp parallel for private(pwrg, beta, c1, c2, c3, c4, c1idx, qidx, c2idx, c3idx, c4idx, betaidx, pwrgidx, jmpterm) firstprivate(J, tempvar, delta, nonminterm, vold, u_opt, u_opt_idx, min_idx)
    for (int rdx = 0; rdx < N_sims; rdx++)
    {
        // This openmp loop does the actual computations and makes use of a linear index from 0 to N_sims-1.
        // For each index in the linear openmp iteration compute
        // the index in the associated grid in order to compute the optimal policy for a specific set of parameters
        qidx = floor(rdx / (pwrggridsteps * c4gridsteps * c3gridsteps * betagridsteps * c2gridsteps * c1gridsteps) % qdists);
        pwrgidx = floor(rdx / (c4gridsteps * c3gridsteps * betagridsteps * c2gridsteps * c1gridsteps) % pwrggridsteps);
        c4idx = floor(rdx / (c3gridsteps * betagridsteps * c2gridsteps * c1gridsteps) % c4gridsteps); // the third inner most
        c3idx = floor(rdx / (betagridsteps * c2gridsteps * c1gridsteps) % c3gridsteps);
        betaidx = floor(rdx / (c2gridsteps * c1gridsteps) % betagridsteps);
        c1idx = floor(rdx / c2gridsteps % c1gridsteps);
        c2idx = rdx % c2gridsteps;

        // Assign to the private variable of the openmp loop the appropriate value of the parameters obtaiend
        // from the index defined above
        c1 = c1grid[c1idx];
        c2 = c2grid[c2idx];
        beta = betagrid[betaidx];
        c3 = c3grid[c3idx];
        c4 = c4grid[c4idx];
        pwrg = pwrggrid[pwrgidx];

        // Here the acatual code for the value iteration alroithm
        J[0] = 0;
        cout << "Starting computation: " << rdx << " out of " << N_sims << "\n";
        do
        {
            // this is do while loop that iterates until the value function stabilises the threshold is delta
            delta = 0;
            for (int i = 1; i < N; i++)
            {
                vold = J[i];
                min_idx = 0;
                for (int udx = 0; udx < noa; udx++)
                { // We check over all possible control the optimal costs and optimal control
                  // depending upon the conidered i the Bellman equation can take different values

                    if (i == N_i)
                    {
                        // in this case, i.e. i = N we only compute the cost as no action is possible to mitigate the spread of the disease since there is a fully infected population
                        tempvar[udx] = (c1 + gamma * i * J[i - 1] + (invdt - gamma * i) * J[i]);
                    }
                    else
                    {
                        // in this case i != N and we compute the cost and the optimal action
                        // note that in the case of burts then two different cases have to considered
                        if (i <= (N_i - qsupportsize))
                        {
                            // a first case is when the support of the burst does not go oever the end of the grid
                            // in the following loop we compute the jump terms for the bellman equation
                            jmpterm = 0;
                            for (int ydx = 0; ydx < qsupportsize; ydx++)
                            {
                                jmpterm += Pqs[qidx][ydx] * J[i + ydx + 1];
                            }
                            // we now add the jump terms to the bellman equation
                            tempvar[udx] = (c1 + c2 * i * pow(u[udx], pwrc2u) + gamma * i * J[i - 1] + c3 * i + c4 * u[udx] + (invdt - gamma * i) * J[i]) + (beta / N_i) * i * (N_i - i) * (pow(1 - 1 * u[udx], pwrg)) * (jmpterm - J[i]);
                        }
                        else
                        { // in this case the support such that the resultin infections
                            // would yield values out of the grid
                            tempvar[udx] = (c1 + c2 * i * pow(u[udx], pwrc2u) + c3 * i + c4 * u[udx] + gamma * i * J[i - 1] + (invdt - gamma * i) * J[i]) + (beta / N_i) * i * (N_i - i) * (pow(1 - 1 * u[udx], pwrg)) * (J[i + 1] - J[i]);
                        }
                    }
                    // in the next if statement the optimal policy is save only if there improvement
                    if (tempvar[udx] < tempvar[min_idx])
                    {
                        min_idx = udx;
                    }
                }
                // save the optimal costs and the optimal action
                J[i] = dt * (*min_element(tempvar, tempvar + noa));
                u_opt[i] = u[min_idx];
                u_opt_idx[i] = min_idx + 1;
                delta = max(delta, abs(vold - J[i]));
            }
        } while (delta > eps);

        // Once the value iteration ends save the data to to two arrays that are then appropriately save to data strcture
        for (int j = 0; j < N; j++)
        {
            ResultsJ[rdx][j] = J[j];
            Resultsu[rdx][j] = u_opt[j];
            Resultsuidx[rdx][j] = u_opt_idx[j];
        }
    }

    // For loop to save the data to the output streams data and thus to file
    cout << "Printing results to file ...";
    for (int ldx = 0; ldx < N_sims; ldx++)
    {
        cout << "\t";
        for (int odx = 0; odx < N; odx++)
        {
            myfile << ResultsJ[ldx][odx] << "\t";
            myfileU << Resultsu[ldx][odx] << "\t";
            myfileUidx << Resultsuidx[ldx][odx] << "\t";
            // cout<< ResultsJ[ldx][odx] << "\n";
            // cout<< Resultsu[ldx][odx] << "\n";
        };
        myfile << std::endl;
        myfileU << std::endl;
        myfileUidx << std::endl;
    }

    myfile.close();
    myfileU.close();
    myfileUidx.close();

    // Now we print to file a summary of the resulting data that is helpful
    // for data analysis/processing
    for (int p = 0; p < c4gridsteps; p++)
    {
        myfileD << c4grid[p] << " ";
    };
    myfileD << std::endl;

    for (int p = 0; p < c3gridsteps; p++)
    {
        myfileD << c3grid[p] << " ";
    };
    myfileD << std::endl;
    
    for (int p = 0; p < betagridsteps; p++)
    {
        myfileD << betagrid[p] / N_i << " ";
    };

    myfileD << std::endl;
    for (int p = 0; p < c1gridsteps; p++)
    {
        myfileD << c1grid[p] << " ";
    };

    myfileD << std::endl;
    for (int p = 0; p < c2gridsteps; p++)
    {
        myfileD << c2grid[p] << " ";
    };
    myfileD << std::endl;

    myfileD << qdists << std::endl;

    for (int p = 0; p < pwrggridsteps; p++)
    {
        myfileD << pwrggrid[p] << " ";
    };
    myfileD << std::endl;

    myfileD << N_i << std::endl;

    for (int p = 0; p < noa; p++)
    {
        myfileD << u[p] << " ";
    };
    myfileD << std::endl;
    myfileD.close();
    return 0;
}
