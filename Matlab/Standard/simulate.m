function [times, values] = simulate(initial_conditions,tspan,rho)
%SIMULATE Simulates the ode with given inital conditions over a time span.
%   Simulate numerically integrates the system of ODEs defined as:
%   R'  = F(R, A)
%   Pwt = G(R, A)
%   A'  = G(A)

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
    
    Aeff = x(3)/x(1);
    dydt = zeros(size(x));
    
    % The change in the ratio of the populations
    if -rho*x(3) + rmax > 0
        dydt(1) = x(1)*(rho*x(3) - rmax*((Aeff^B1)/(K1^B1 + Aeff^B1)) - gammawt*Aeff/(K1 + Aeff));
    else
        dydt(1) = x(1)*(rmax*(1 - (Aeff^B1)/(K1^B1 + Aeff^B1)) - gammawt*Aeff/(K1 + Aeff));
    end
    
    % The change in the proportion of the WT present now compared to the
    % start
    dydt(2) = rmax*(1 - (Aeff^B1)/(K1^B1 + Aeff^B1)) - gammawt*Aeff/(K1 + Aeff)*x(2);
    
    % Change in the antibiotic concentration
    if x(4) == 1
        dydt(3) = alpha - gamma*x(3);
    else
        dydt(3) = -gamma*x(3);
    end
end
end

