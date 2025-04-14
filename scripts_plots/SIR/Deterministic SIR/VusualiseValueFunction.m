
% This function displays value function (for SIR model) for the the S and I compartmetn at time step t0

filename = 'valueFunction.t0';
fid = fopen(filename,'r');

% Read the first line of the file(ModeNB, State0, State1...)
C = textscan(fid,'%u %u %u %u %u',1,'Delimiter','\n');

% Identify differents parameters
ModeNB = C{1};
StateDim = C{2};

for i = 1:StateDim
    X(i) = C{2+i}; % get the state dimensions for the two variables
end

A = importdata(filename,'\t',1); 
v = A.data;

for j = 1:ModeNB
    switch StateDim
        case 1
            disp('No need reshape');
            figure;
            plot(v(:,j));
        case 2
            V1 = reshape(v(:,j), X(1), X(2));
            V1 =V1.*(V1<1e10); % this command removes the the high values in the optimal value funciton 
            % that are generated on the boundary
            figure;
            p=surf(V1);
            set(p,'edgecolor','none');
            xlabel('I');
            ylabel('S');
            c = colorbar;
            colormap (jet(1024));
            ylabel(c,'Value function');
        otherwise
            warning('Error! More than 2 dimensions')
    end
end

fclose(fid);