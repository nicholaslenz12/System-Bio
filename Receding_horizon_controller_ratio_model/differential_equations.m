function values = DifferentialEquations(initialConditions,tspan,rho)
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
% values            : The value of both popoulation ratios and the
%                     antibiotic concentration at each time step.
%% Important Parameters

rmax = log(2)/20;
B1 = 1.7;
K1 = 10;
gamma = log(2)/20;
gammawt = .1;
alpha = 1;
%% Hill Function

function hill_output = hill_function(half_occupation,hill_coefficient,concentration)
    hill_output = (concentration^hill_coefficient) ...
                 /(half_occupation^hill_coefficient + concentration^hill_coefficient);
end
%% Time Derivative

function dvdt = growth_control(t, x)
    % Sets the effective antibiotic concentration, which is Pc*A/Pwt or
    % A/R. The vector of the derivatives is pre-allocated.
    Aeff = x(3)/x(1);
    dvdt = zeros(size(x));
    
    % Computes R'
    if -rho*x(3) + rmax > 0
        dvdt(1) = x(1)*(rho*x(3) - rmax*hill_function(K1,B1,Aeff) - gammawt*hill_function(K1,1,Aeff));
    else
        dvdt(1) = x(1)*(rmax*(1 - hill_function(K1,B1,Aeff)) - gammawt*hill_function(K1,1,Aeff));
    end
    
    % Computes P'wt
    dvdt(2) = (rmax*(1 - hill_function(K1,B1,Aeff)) - gammawt*hill_function(K1,1,Aeff))*x(2);
    
    % Computes A'
    if x(4) == 1
        dvdt(3) = alpha - gamma*x(3);
    else
        dvdt(3) = -gamma*x(3);
    end
end
%% Simulate ODE

[~, values]=ode45(@growth_control,tspan,initialConditions);

end

