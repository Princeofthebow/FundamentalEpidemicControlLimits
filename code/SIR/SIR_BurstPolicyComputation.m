clear all
clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The code in this file computes the optimal policies using value iteration 
% for the Susceptible-Infected-Recovered with bursts model via the value iteration algoritihm for 
% the considered costs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N_i = 50;% Number of possible infective
N = N_i + 1;% number of nodes in markov chain
eps = 1e-6;% Threshold for convergence of the Value iteration algorithm

nu = 1e2;% uniformisation factor denoted as nu in the paper
% The uniformisation factor is here selected empirically to be greater than
% the sum of the rates exiting a node. The article quantifies this
% precisely but we use empirical values in the code that have to be adapted
% accordinlgy based upon N

% System constants
gamma = 0.32;
R0 = 3.5;
mu = R0 * gamma / N_i;
ub = linspace(0, 0.8, 15);

% Definition of cost parameters as in the paper.
% The costs function is c1+c2ih(u)+c3i+c4Nu
c1 = 0;
c2 = 1;
c3 = 1;
c4 = 0;
qeff = .8;
pwrh = 3;
% The computation also consider h(u) = (1-qeff*u)^pwrh and simple adaptations are
% required to consider other h(u)


% Burst infection distribution
Nq = 3; % maximum burst size (Nb has to be smaller than N)
Q = ones(1, Nq) / Nq;  % uniform burst distribution

% Initalisation of the data structures for the value iteration algorithm
minpos = NaN(N);
minpos(:,1) = ones(N,1);
J = zeros(N, N);% array holding the value function 
J(:,1) = 0 * 1/nu * ones(N,1);  % Initial values for J

tic
while true
    delta = 0;
    for sdx = 1:N_i 
        s = sdx - 1;
        for idx = 2:N - s  
            i = idx - 1;
            Vnew = J(sdx, idx);
            switch s
                case 0
                    [minval, minpos(sdx, idx)] = min(c1 + c2 * i * ub + c3 * i + c4 * ub.*N);
                    J(sdx, idx) = (1 / nu) * (gamma * i * J(sdx, idx - 1) + (nu - gamma * i) * J(sdx, idx) + minval);
                otherwise
                    costs = zeros(size(ub));
                    for uidx = 1:length(ub)
                        u = ub(uidx);
                        sumQ = 0;
                        for k = 1:Nq
                            if s - k >= 0 && i + k <= N_i
                                s_next = s - k;
                                i_next = i + k;
                                sumQ = sumQ + Q(k) * J(s_next + 1, i_next + 1);
                            else
                                sumQ = sumQ + Q(k) * J(sdx, idx); 
                            end
                        end
                        transition_term = mu * (1 - qeff * u)^pwrh * s * i * (sumQ - J(sdx, idx));
                        costs(uidx) = c1 + c2 * i * u + c3 * i + c4 * u*N + transition_term;
                    end
                    [minval, minpos(sdx, idx)] = min(costs);
                    J(sdx, idx) = (1 / nu) * (gamma * i * J(sdx, idx - 1) + (nu - gamma * i) * J(sdx, idx) + minval);
            end
            delta = max(delta, abs(Vnew - J(sdx, idx)));
        end
    end
delta
    if delta < eps
        break;
    end
end
etime = toc;

disp(['The computation took '  num2str(etime)  ' seconds'])

figure();
surf(0:N_i,0:N_i,J)
xlabel('i (infected)')
ylabel('s (susceptible)')
zlabel('J(s,i)')
title('Value Function for the SIR model with bursts')

policy = minpos;
policy(~isnan(policy)) = ub(policy(~isnan(policy)));

figure();
surf(0:N_i, 0:N_i, policy)
xlabel('I (infected)')
ylabel('S (susceptible)')
zlabel('Optimal control')
title('Optimal Policy with Burst Infections')