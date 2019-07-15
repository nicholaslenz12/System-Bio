function [values] = model_2(initialConditions,tspan)
%% MODEL_2 Simulates the ODE:
%  Pwt' = f(Pwt, Pc, A)
%  Pc'  = g(Pwt, Pc, A)
%  A'   = h(A)
%Args
%  initialConditions
%      Initial values for Pwt, Pc, A
%  tspan
%      Length of simulation time.
%Rets
%  values
%      The state of Pwt, Pc, A at given time steps in the simulation

%% ------------------------------------------------------------------------
% PARAMETERS
% -------------------------------------------------------------------------
r1 = log(2)/20;
r2 = log(2)/40;
alpha = 0.01;
K1 = 0.01;
K2 = 3;
B1 = 1.5;
B2 = 1.5;
gammawt = 0.01*r1;
gammac = 0.001*r2;
gamma = log(2)/40;
S = 30;

function saturation = hill_function(half_occupation,hill_coefficient,concentration)
%% HILL_FUNCTION
%  Computes the value of the hill function for specific parameters
%Args
%  half_occupation
%      The half occupation. In otherwords, when the
%      concentration of molecule reaches half_occupation, the
%      hill output is exactly 0.5.
%  hill_coefficient
%      cooperative binding coefficient.
%  concentration
%      concentration of the molecule.
%Rets
%  saturation
%      saturation between [0,1].
    saturation = (concentration^hill_coefficient) ...
                 /(half_occupation^hill_coefficient + concentration^hill_coefficient);
end

%% ------------------------------------------------------------------------
% DERIVATIVE
% -------------------------------------------------------------------------
function dxdt = growth_control(t, x)
%% GROWTH_CONTROL
%  Computes the Pwt', Pc', and A'
%Args
%  x
%      vector of state variables x = (Pwt, Pc, A, u)
%Rets
%
    Pwt = x(1);
    Pc  = x(2);
    A   = x(3);
    u   = x(4);
    Aeff = Pc*A/Pwt; % Effective concentration of the antibiotic.
    dxdt = zeros(size(x)); % Initializes the derivative.

    % Computes P'wt
    dxdt(1) = r1*(1 - hill_function(K1, B1, Aeff))*Pwt*(1 - (Pwt+Pc)/(Pwt+Pc + S)) ...
              - gammawt*Aeff*Pwt;

    % Computes P'c
    dxdt(2) = r2*(1 - hill_function(K2, B2, A))*Pc*(1 - (Pwt+Pc)/(Pwt+Pc + S)) ...
              - gammac*Pc;

    % Computes A for chosen proportions of alpha.
    if x(4) == 2
        dxdt(3) = alpha/2 - gamma*x(3);
    elseif x(4) == 3
        dxdt(3) = alpha - gamma*x(3);
    else
        dxdt(3) = -gamma*x(3);
    end
end

%% ------------------------------------------------------------------------
% MAIN
% -------------------------------------------------------------------------
[~, values]=ode45(@growth_control,tspan,initialConditions);

end
