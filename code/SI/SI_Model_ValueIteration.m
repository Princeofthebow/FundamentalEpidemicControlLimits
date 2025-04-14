clear all
clc

N_i = 20; % Number of possible infective
N = N_i+1; % number of nodes in markov chain
epsilon = 1e-6; % Threshold for convergence of the Value iteration algorithm

nu = 100; % uniformisation factor denoted as nu in the paper
dt = 1/nu; 

% System constants
gamma = .32; % Recovey rate constant
beta = 1.9*gamma; % Infection rate constant
u = linspace(0,0.8,2); % set of available control action

% Definition of cost parameters
c1 = 1;
c2 = 1;
c3 = 0;
c4 = 0;

% Initalisation of the data sctructuers for the value iteration algorithm
J = zeros(N,1);  % array holding the value iteration                  
%J(1) = dt*c1;
minpos(1) = length(u);

tic
while true

    delta = 0;

    for i = 1:N_i
        idx = i+1;                   

        Vold = J(idx);
        switch i
            case N_i
                   [minval(idx), minpos(idx)] = min(c1);
                   J(idx) = dt*(gamma*i.*J(idx-1) + (1/dt-gamma*i).*J(idx) );              
            otherwise
                   [minval(idx), minpos(idx)] = min( (c1 + c2.*i.*u + c3*i + c4*u) + beta./N_i.*(1-u).*i.*(N_i-i).*( J(idx+1)-J(idx) ) );
                   J(idx) = dt*(gamma*i.*J(idx-1) + (1/dt-gamma*i).*J(idx)+ minval(idx)  );
        end 

        delta = max(delta,abs(Vold-J(idx)));      
    end

    if delta < epsilon
        break;
    end

    
end
etime = toc
J
