function values = differential_equations(initialConditions,tspan)
%DIFFERENTIALEQUATIONS Simulates the ode:
%
% over a given timespan tspan.
%The inputs for this function are:
% initialConditions : The initial state of the system.
% tspan             : A vector of time steps over which the simulation
%                     runs.
%The outputs are:
% values            : The value of both popoulation ratios and the
%                     antibiotic concentration at each time step.
%% Important Parameters ---------------------------------------------------

rmax = log(2)/20;
B1 = 1.7;
K1 = 5;
K2 = 1;
gamma = .1;
gammawt = .1;
alpha = 70;
S = 60;
%% Hill Function ----------------------------------------------------------
%HILL FUNCTION computes the value of the hill function for specific
% parameters
%The inputs for this function are:
% half_occupation  : The half occupation. In otherwords, when the
%                  concentration of molecule reaches half_occupation, the
%                  hill output is exactly 0.5.
% hill_coefficient : cooperative binding coefficient.  
% concentration    : concentration of the molecule.
%The outputs are:
% hill_output      : saturation between [0,1].

function hill_output = hill_function(half_occupation,hill_coefficient,concentration)
    hill_output = (concentration^hill_coefficient) ...
                 /(half_occupation^hill_coefficient + concentration^hill_coefficient);
end
%% Time Derivative --------------------------------------------------------

function dxdt = growth_control(t, x)
%GROWTH CONTROL Computes the time derivative of P'wt, P'c, and A for state x.
%The inputs for this function are:
% t : unused!
% x : The state of the system where:
%     x(1) == P'wt
%     x(2) == P'c
%     x(3) == A

    Aeff = x(2)*x(3)/x(1); % Computes the effective concentration of the antibiotic.
    
    dxdt = zeros(size(x)); % Initializes the time derivative vector.
    
    % Computes P'wt
    if x(1) > .0001
        dxdt(1) = (rmax*(1-hill_function(K1,B1,Aeff)))*x(1)*(1-(x(1)+x(2))/(x(1)+x(2)+ S))-gammawt*Aeff;
    else
        dxdt(1) = 0;
    end
        
    % Computes P'c
    if x(2) > .0001
        dxdt(2) = (rmax*(1 - hill_function(K2,B1,x(3))))*x(2)*(1-(x(1)+x(2))/(x(1)+x(2)+ S));
    else
        dxdt(1) = 0;
    end
        
    % Computes A for chosen proportions of \alpha.
    dxdt(3) = alpha - gamma*x(3);
end
%% Simulate ODE

[~, values]=ode45(@growth_control,tspan,initialConditions);

end