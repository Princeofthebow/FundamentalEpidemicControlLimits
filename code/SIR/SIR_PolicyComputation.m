clear all
clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The code in this file computes the optimal policies using value iteration 
% for the Susceptible-Infected-Recovered model via the value iteration algoritihm for 
% the considered costs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N_i = 15;% Number of possible infective
N = N_i+1; % number of nodes in markov chain
eps = 1e-6;% Threshold for convergence of the Value iteration algorithm

nu = 1e2;% uniformisation factor denoted as nu in the paper
% The uniformisation factor is here selected empirically to be greater than
% the sum of the rates exiting a node. The article quantifies this
% precisely but we use empirical values in the code that have to be adapted
% accordinlgy based upon N

% System constants
gamma = 0.32;
R0= 3.5;
mu = R0*gamma/N_i;
ub = linspace(0,0.8,2);

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
minpos = NaN(N);
minpos(:,1) = ones(N_i+1,1);
J = zeros(N_i+1,N_i+1);% array holding the value function  
J(:,1)=0*1/nu*ones(N_i+1,1);

tic
while true
    
    delta = 0;
    for sdx = 1:N_i
    s = sdx-1;
        for idx = 2:N_i+1-s
            i=idx-1;
            Vnew = J(sdx,idx);
            switch s
                    case 0
                         [minval,minpos(sdx,idx)] = min( c1+c2.*i.*ub+c3.*i+c4.*ub.*N);
                        J(sdx,idx) = (1./nu).*( gamma*i.*J(sdx,idx-1) + (nu-gamma*i).*J(sdx,idx) + minval);  % as function?

                    otherwise
                        [minval,minpos(sdx,idx)] = min( c1+c2.*i.*ub+c3.*i+c4.*ub.*N+ mu.*(1-ub).*s.*i.*( J(sdx-1,idx+1)-J(sdx,idx) ) );
                        J(sdx,idx) = (1./nu).*( gamma*i.*J(sdx,idx-1) + (nu-gamma*i).*J(sdx,idx) + minval);  % as function?
            end            
            delta = max(delta,abs(Vnew-J(sdx,idx)));  
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
title('Value Function for the SIR model')


policy = minpos;
policy(find(~isnan(policy))) = ub((minpos(find(~isnan(policy)))));

figure();
surf(0:N_i,0:N_i,policy)
xlabel('i (infected)')
ylabel('s (susceptible)')
zlabel('u^*(s,i)')
title('Optimal policy for the SIR model')
