clear all
clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The code in file runs the Stochastic simualtion algorithm for the
%  Susceptible-Infected model.
% The code does not implement any particular control policy; minor edits
% are required for this
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

gamma= 0.32; % reovery rate
N = 10; % Number of infected individuals
betaf = 1.6*gamma/N;% infection rate constant

iSimTime = 15; % Simulation time 

%%%%%
iReal = 2; % number of realisation
Data = {}; % data structure for results
Time = {}; % data structure hodling the time of the events

% State change matrix
WMatrix = [+1, -1];

tic  
parfor idx = 1:iReal
    t=0;
    x= [5];
    xf=x;
    tf = t;

    disp(['Working on realisation: ' num2str(idx)])    

    while t<=iSimTime
             
       rates = [ betaf*x(1)*(N-x(1)), gamma*x(1)];  

       % a particular policy can be implemented by changing the rate a this
       % point 

        ratesum = sum(rates);
        
        tau=-log(rand)/ratesum; % can be changed to ;
        
        reac = sum(cumsum(rates/ratesum)<rand)+1;

        x = x + WMatrix(:,reac);

        t = t + tau;
        tf = [tf, t];
        xf = [xf,x];

    end
    Data{idx} = xf;
    Time{idx} = tf;
end

etime = toc;
disp(['The simulation took: ' num2str(etime) ' seconds' ])  

figure();
for idx = 1:iReal
hold on
stairs(Time{idx},Data{idx});
end
ylabel('i (infected)')
xlabel('y')
title('Time realisation of the SI model')