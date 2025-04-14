clear all
clc

N_i = 50;
N = N_i+1;
eps = 1e-6;
c1= 0;
c2= 10; 
c3 = 1;

nu = 1e2;

gammai = 0.32;
gammae= 0.5;
mu = 3.5*gammai/N_i;
ub = linspace(0,0.8,10);

J=zeros(N,N,N);
minpos = NaN(N,N,N);


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

figure()
stoplot = [3,4];%floor(linspace(2,N-30,4));
for udx = 1:length(stoplot)
subplot(length(stoplot),1,udx)
% surf(0:N_i,0:N_i,ub(squeeze(minpos(:,stoplot(udx),:))))
surf(0:N_i,0:N_i,squeeze(policy(:,stoplot(udx),:)))
%effectR0 = 1.6*(1-DataU{tdx}');
view(0,90)
end
