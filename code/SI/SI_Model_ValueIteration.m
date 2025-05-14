clear all
clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The code in this file computes the optimal policies using value iteration 
% for the Susceptible-Infected model via the value iteration algoritihm for 
% the considered costs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N_i = 20; % Number of possible infective
N = N_i+1; % number of nodes in markov chain
epsilon = 1e-6; % Threshold for convergence of the Value iteration algorithm

nu = 100; % uniformisation factor denoted as nu in the paper
% The uniformisation factor is here selected empirically to be greater than
% the sum of the rates exiting a node. The article quantifies this
% precisely but we use empirical values in the code that have to be adapted
% accordinlgy based upon N

% System constants
gamma = .32; % Recovey rate constant
beta = 1.9*gamma; % Infection rate constant
u = linspace(0,0.8,2); % set of available control action

% Definition of cost parameters as in the paper.
% The costs function we consider is c1+c2iz(u)+c3i+c4Nu
% The computation consider z(u) = u and simple adaptations are
% required to consider other z(u)
c1 = 0;
c2 = 1;
c3 = 1;
c4 = 0;
% In this computation we consider h(u) = 1-u but with minor edits other
% functions can be considered


% Initalisation of the data structures for the value iteration algorithm
J = zeros(N,1);  % array holding the value function                 
minpos(1) = length(u); 
% arbitrary value as with no infected there is not mitigation measure to be taken

tic
while true

    delta = 0;

    for i = 1:N_i
        idx = i+1;                   

        Vold = J(idx);
        switch i
            case N_i
                   [minval(idx), minpos(idx)] = min(c1);
                   J(idx) = 1/nu*(gamma*i.*J(idx-1) + (nu-gamma*i).*J(idx) );              
            otherwise
                   [minval(idx), minpos(idx)] = min( (c1 + c2.*i.*u + c3*i + c4*N*u) + beta./N_i.*(1-u).*i.*(N_i-i).*( J(idx+1)-J(idx) ) );
                   J(idx) = 1/nu*(gamma*i.*J(idx-1) + (nu-gamma*i).*J(idx)+ minval(idx)  );
        end 

        delta = max(delta,abs(Vold-J(idx)));      
    end

    if delta < epsilon
        break;
    end

    
end

etime = toc;
disp(['The computation took '  num2str(etime)  ' seconds'])
figure();
stairs(0:N_i, J)
xlabel('i (infected)')
ylabel('J(i)')
title('Value Function for the SI model')

figure();
stairs(1:N_i-1, u(minpos(2:end-1)))
xlabel('i (infected)')
ylabel('u^*(i)')
title('Optimal control for the SI model')
