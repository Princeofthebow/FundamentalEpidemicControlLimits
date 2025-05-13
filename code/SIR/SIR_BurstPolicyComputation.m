clear all
clc

N_i = 80;
N = N_i + 1;
eps = 1e-6;
J = zeros(N, N);
c1 = 0;
c2 = 1;
c3 = 1;
c4 = 0;
nu = 1e2;
c2pwr = 2;
c4pwr = 0;
c3pwr = 1;
pwrh = 3;
qeff = .8;

gamma = 0.32;
R0 = 3.5;
mu = R0 * gamma / N_i;
ub = linspace(0, 0.8, 15);

% Burst infection distribution
Nq = 3; % maximum burst size
Q = ones(1, Nq) / Nq;  % uniform burst distribution

minpos = NaN(N);
minpos(:,1) = ones(N,1);
J(:,1) = 0 * 1/nu * ones(N,1);  % Initial values for J

while true
    delta = 0;
    for sdx = 1:N_i  % s = 0 to N_i - 1
        s = sdx - 1;
        for idx = 2:N - s  % i = 1 to N_i - s
            i = idx - 1;
            Vnew = J(sdx, idx);
            switch s
                case 0
                    [minval, minpos(sdx, idx)] = min(c1 + c2 * i * ub.^c2pwr + c3 * i^c3pwr + c4 * ub.^c4pwr);
                    J(sdx, idx) = (1 / nu) * (gamma * i * J(sdx, idx - 1) + (nu - gamma * i) * J(sdx, idx) + minval);
                otherwise
                    costs = zeros(size(ub));
                    for uidx = 1:length(ub)
                        u = ub(uidx);
                        sumQ = 0;
                        for k = 1:Nq
                            if s - k >= 0 && i + k <= N_i
                                s_next = s - k;
                                i_next = i + k;
                                sumQ = sumQ + Q(k) * J(s_next + 1, i_next + 1);
                            else
                                sumQ = sumQ + Q(k) * J(sdx, idx);  % fallback: stay in current
                            end
                        end
                        transition_term = mu * (1 - qeff * u)^pwrh * s * i * (sumQ - J(sdx, idx));
                        costs(uidx) = c1 + c2 * i * u^c2pwr + c3 * i^c3pwr + c4 * u^c4pwr + transition_term;
                    end
                    [minval, minpos(sdx, idx)] = min(costs);
                    J(sdx, idx) = (1 / nu) * (gamma * i * J(sdx, idx - 1) + (nu - gamma * i) * J(sdx, idx) + minval);
            end
            delta = max(delta, abs(Vnew - J(sdx, idx)));
        end
    end
delta
    if delta < eps
        break;
    end
end

figure();
surf(0:N_i, 0:N_i, J)
xlabel('I (infected)')
ylabel('S (susceptible)')
zlabel('Cost-to-go')
title('Value Function with Burst Infections')

policy = minpos;
policy(~isnan(policy)) = ub(policy(~isnan(policy)));

figure();
surf(0:N_i, 0:N_i, policy)
xlabel('I (infected)')
ylabel('S (susceptible)')
zlabel('Optimal control')
title('Optimal Policy with Burst Infections')
