clear all
clc

N_i = 15;
N = N_i+1;
eps = 1e-6;
J = zeros(N_i+1,N_i+1);
c1 = 0;
c2 = 1;
c3 = 1;
c4 = 0;
nu = 1e2;
c2pwr = 1;
c4pwr = 1;
c3pwr = 1;
pwrh = 1;
qeff = 1;

gamma = 0.32;
R0= 3.5;
mu = R0*gamma/N_i;
ub = linspace(0,0.8,15);

minpos = NaN(N);
minpos(:,1) = ones(N_i+1,1);
J(:,1)=0*1/nu*ones(N_i+1,1);
while true
    
    delta = 0;
    for sdx = 1:N_i %si 1:N+1, s = 0:N
    s = sdx-1;
        for idx = 2:N_i+1-s
            i=idx-1;
            Vnew = J(sdx,idx);
            switch s
                    case 0
                         [minval,minpos(sdx,idx)] = min( c1+c2.*i.*ub.^c2pwr+c3.*i.^c3pwr+c4.*ub.^c4pwr);
                        J(sdx,idx) = (1./nu).*( gamma*i.*J(sdx,idx-1) + (nu-gamma*i).*J(sdx,idx) + minval);  % as function?

                    otherwise
                        [minval,minpos(sdx,idx)] = min( c1+c2.*i.*ub.^c2pwr+c3.*i.^c3pwr+c4.*ub.^c4pwr+ mu.*(1-qeff*ub).^pwrh.*s.*i.*( J(sdx-1,idx+1)-J(sdx,idx) ) );
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
figure();
surf(0:N_i,0:N_i,J)

policy = minpos;
policy(find(~isnan(policy))) = ub((minpos(find(~isnan(policy)))));

figure();
surf(0:N_i,0:N_i,policy)

% figure();
% Re=R0*linspace(0,1,N)'.*ones(N).*(1-policy);
% surf(0:N_i,0:N_i,Re)
% HIT= (1-1/R0)*N_i
% 
% ExportSurfplot('SIR_MJP_N100c2c3_10.dat',0:N_i,0:N_i,policy',1,1)