clear all
clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The code in file runs the Stochastic simualtion algorithm for the
%  Susceptible-Infected-Recovered model.
% The code does not implement any particular control policy; minor edits
% are required for this
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



gamma = 0.32; % Recovery rate
N_i = 500;%  Number of individuals in the population
N = N_i+1;% individuals in the grid
R0 = 2.5; % R0 of the disease
mu = R0*gamma/N_i; % this is the total infection rate constant 

iSimTime = 60; % Maximum Simulation Time
iReal = 5; % number of realisation to be considered
Data = {}; % Data structure for the trajectories
Time={}; % Data structure for the trajectories

% state change matrix 
WMatrix = [-1, 0; 
           1,-1];

init_infec = 3;% initial level of infection

tic % tic is for starting a clock to obtain the computation time
parfor idx = 1:iReal
    % this loop is for realisations

    t=0;
    x= [N_i-init_infec;init_infec];

    xf=x;
    tf = t;

    disp(['Working on realisation: ' num2str(idx)])    

    while t<=iSimTime % from here on this is standard Gillespie algorithm
             
       rates = [ mu*x(1)*x(2), gamma*x(2)]; % here we compute the rates
       
       % a particular policy can be implemented by changing the rate a this
       % point 
       %eg
       % if (x(2)>50)
       % rates = [ 0.65*mu*x(1)*x(2), gamma*x(2)];
       % end

        ratesum = sum(rates); % sum of rates
        
        tau=-log(rand)/ratesum; % extraction of next time
        
        reac = sum(cumsum(rates/ratesum)<rand)+1; % computation fo next reaction

        x = x + WMatrix(:,reac); % state change

        t = t + tau;
        tf = [tf, t];
        xf = [xf,x];

    end
    % save data to output structures
    Data{idx} = xf;
    Time{idx} = tf;

end

etime = toc; % stop time for execution time
disp(['It took: ' num2str(etime) ' seconds' ])  
etime = toc;
disp(['The simulation took: ' num2str(etime) ' seconds' ])  

figure();
for idx = 1:iReal
plot3(Data{1,idx}(1,:),(Time{idx}),(Data{1,idx}(2,:)),'b')
hold on
end
zlabel('i (infected)')
xlabel('s (susceptible)')
ylabel('t')
title('Time realisation of the SIR model')