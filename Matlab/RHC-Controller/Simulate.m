function [times, values] = Simulate(initialConditions,tspan,rho)
%SIMULATE Simulates the ode:
%   R'  = F(R, A)
%   P'wt = G(R, A)
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
%% Important Parameters

rmax = log(2)/20;
B1 = 1.7;
K1 = 10;
gamma = log(2)/20;
gammawt = .1;
alpha = 1;
%% Simulate ODE

[times, values]=ode45(@growth_control,tspan,initialConditions);

function dvdt = growth_control(t, x)
    % Sets the effective antibiotic concentration, which is Pc*A/Pwt or
    % A/R. The vector of the derivatives is pre-allocated.
    Aeff = x(3)/x(1);
    dvdt = zeros(size(x));
    
    % Computes R'
    if -rho*x(3) + rmax > 0
        dvdt(1) = x(1)*(rho*x(3) - rmax*((Aeff^B1)/(K1^B1 + Aeff^B1)) - gammawt*Aeff/(K1 + Aeff));
    else
        dvdt(1) = x(1)*(rmax*(1 - (Aeff^B1)/(K1^B1 + Aeff^B1)) - gammawt*Aeff/(K1 + Aeff));
    end
    
    % Computes P'wt
    dvdt(2) = rmax*(1 - (Aeff^B1)/(K1^B1 + Aeff^B1)) - gammawt*Aeff/(K1 + Aeff)*x(2);
    
    % Computes A'
    if x(4) == 1
        dvdt(3) = alpha - gamma*x(3);
    else
        dvdt(3) = -gamma*x(3);
    end
end
end

