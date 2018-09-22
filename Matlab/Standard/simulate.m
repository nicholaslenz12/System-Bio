function [times, values] = simulate(initialConditions,tspan,rho)
%SIMULATE Simulates the ode:
%   R'  = F(R, A)
%   Pwt = G(R, A)
%   A'  = G(A)
% over a given timespan tspan, with metabolic cost rho.
%The inputs for this function are:
% initialConditions : The initial state of the system.
% tspan             : A vector of time steps over which the simulation
%                     runs.
% rho               : The metablic cost on the controller for producing
%                     antibiotic.
%The outputs are:
% times             : A vector of time steps over which the simulation
%                     runs. < Why do you have this?
% values            : The value of both popoulation ratios and the
%                     antibiotic concentration at each time step.
rmax = log(2)/20;
B1 = 1.7;
K1 = 10;
gamma = log(2)/20;
gammawt = .1;
alpha = 1;

%% Simulate ODE
[times, values]=ode45(@growth_control,tspan,initialConditions);

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

