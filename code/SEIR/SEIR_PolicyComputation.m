clear all
clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The code in this file computes the optimal policies using value iteration 
% for the Susceptible-Exposed-Infected-Recovered model 
% via the value iteration algoritihm for the considered costs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N_i = 50;% Number of possible infective
N = N_i+1;% number of nodes in markov chain
eps = 1e-6;% Threshold for convergence of the Value iteration algorithm

nu = 1e2;% uniformisation factor denoted as nu in the paper
% The uniformisation factor is here selected empirically to be greater than
% the sum of the rates exiting a node. The article quantifies this
% precisely but we use empirical values in the code that have to be adapted
% accordinlgy based upon N

% System constants
gammai = 0.32;
gammae= 0.5;
mu = 3.5*gammai/N_i;
ub = linspace(0,0.8,5);

% Definition of cost parameters as in the paper.
% The costs function is c1+c2iz(u)+c3i
% The computation consider z(u) = u and simple adaptations are
% required to consider other z(u)
c1= 0;
c2= 10; 
c3 = 1;
% In this computation we consider h(u) = 1-u but with minor edits other
% functions can be considered

% Initalisation of the data structures for the value iteration algorithm
J=zeros(N,N,N);% array holding the value function
minpos = NaN(N,N,N);

tic
while true
    
    delta = 0;
    
    for sdx = 1:N_i+1
    s = sdx-1;
        for edx = 2:N_i+1 -s
         e = edx-1;
            for idx = 2:N_i+1-s-e
                i=idx-1;
                Vnew = J(sdx,edx,idx);
                switch s
                        case 0
                             [minval,minpos(sdx,edx,idx)] = min( c1+c2.*i.*ub+ c3.*i);
                            J(sdx,edx,idx) = (1./nu).*( gammai*i.*J(sdx,edx,idx-1) + gammae*e.*J(sdx,edx-1,idx+1)+ (nu-gammai*i-gammae*e).*J(sdx,edx,idx) + minval);  % as function?
    
                        otherwise
                            [minval,minpos(sdx,edx,idx)] = min( c1+c2.*i.*ub+ c3.*i+ mu.*(1-ub).*s.*i.*( J(sdx-1,edx+1,idx)-J(sdx,edx,idx) ) );
                            J(sdx,edx,idx) = (1./nu).*( gammai*i.*J(sdx,edx,idx-1) + gammae*e.*J(sdx,edx-1,idx+1)+ (nu-gammai*i-gammae*e).*J(sdx,edx,idx) + minval);  % as function?
                end
                delta = max(delta,abs(Vnew-J(sdx,edx,idx)));  
            end
        end
    end

    delta
    
    if delta < eps
        break;
    end
   policy = minpos;
policy(find(~isnan(policy))) = ub((minpos(find(~isnan(policy)))));
  
end

etime = toc;
disp(['The computation took '  num2str(etime)  ' seconds'])

figure()
stoplot = [3,4];
for udx = 1:length(stoplot)
subplot(length(stoplot),1,udx)
surf(0:N_i,0:N_i,squeeze(policy(:,stoplot(udx),:)))
view(0,90)
end
title('Optimal policy for the SEIR model')
