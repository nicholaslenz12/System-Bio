function [times, values] = differential_equations(initial_conditions,tspan,rho)
%SIMULATE Simulates the ode with given inital conditions over a time span.
%   Simulate numerically integrates the system of ODEs defined as:
%   WT' = F(WT, C, A)
%   C   = G(A)
%   A   = H(WT, C, A)

% values
rmax = log(2)/20;
B1 = 1.7;
K1 = 10;
gamma = log(2)/20;
gammawt = .1;
alpha = 1;

%% Simulate ODE
[times, values]=ode45(@growth_control,tspan,initial_conditions);

function dydt = growth_control(t, x)
    Aeff = x(2)*x(3)/x(1);
    
    dydt = zeros(size(x));
    dydt(1) = rmax*(1 - (Aeff^B1)/(K1^B1 + Aeff^B1))*x(1) ...
              - (gammawt*Aeff/(K1 + Aeff))*x(1);
    if -rho*x(3) + rmax > 0
        dydt(2) = (-rho*x(3) + rmax)*x(2);
    else
        dydt(2) = 0;
    end
    if x(4) == 1
        dydt(3) = alpha - gamma*x(3);
    else
        dydt(3) = -gamma*x(3);
    end
end
end

