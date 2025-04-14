%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%      File for importing C++ results to Matlab    %%%%%%%%%%%%%%
%%%%%%%%%%%                SI Model                          %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
clc

dataJ = load('SI_VI.txt'); % load the data for the optimal value function
dataU = load('SI_VI_U.txt'); % load the data for the optimal policy
dataD = importdata('RunData.txt',' '); % load the summary for the computations
% This last file defines various variable such as the size of the considered grids as
% well the size of the consideted Markov Chain. See more details below


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Here we import the data into a table then copy it into an array for ease
% of processing.
% The order of the varaibles and parameter grids written to file is
% c4grid
% c3grid
% betagrid
% c1grid
% c2grid
% pwggrid
% ub
% N_i
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ReadFileTable = readtable('RunData.txt');

c4grid = ReadFileTable{1,:}( ~isnan(ReadFileTable{1,:}) );
c4gridsize = length(c4grid);

c3grid = ReadFileTable{2,:}( ~isnan(ReadFileTable{2,:}) );
c3gridsize = length(c3grid);

betagrid = ReadFileTable{3,:}( ~isnan(ReadFileTable{3,:}) );
betagridsize = length(betagrid);

c1grid = ReadFileTable{4,:}( ~isnan(ReadFileTable{4,:}) );
c1gridsize = length(c1grid);

c2grid = ReadFileTable{5,:}( ~isnan(ReadFileTable{5,:}) );
c2gridsize = length(c2grid);


pwggrid = ReadFileTable{6,:}( ~isnan(ReadFileTable{6,:}) );
pwggridsize = length(pwggrid);

N_i = ReadFileTable{7,:}( ~isnan(ReadFileTable{7,:}) );

ub = ReadFileTable{8,:}( ~isnan(ReadFileTable{8,:}) );

NSim = pwggridsize*c4gridsize*c3gridsize*betagridsize*c1gridsize*c2gridsize;

for i=0:NSim-1 % from 0 to Number of Simulations
    pwgidx = floor(mod(i/(c4gridsize*c3gridsize*betagridsize*c1gridsize*c2gridsize),pwggridsize))+1;
    c4idx = floor(mod(i/(c3gridsize*betagridsize*c1gridsize*c2gridsize),c4gridsize))+1;
    c3idx = floor(mod(i/(betagridsize*c1gridsize*c2gridsize),c3gridsize))+1;
    betaidx = floor(mod(i/(c1gridsize*c2gridsize),betagridsize))+1;
    c1idx = floor(mod(i/c2gridsize,c1gridsize))+1;
    c2idx = mod(i,c2gridsize)+1;
    
    J_opt{c1idx,c2idx,betaidx,c3idx,c4idx,pwgidx} = dataJ(i+1,:);
    u_opt{c1idx,c2idx,betaidx,c3idx,c4idx,pwgidx} = dataU(i+1,:);
end