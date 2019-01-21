%MAIN Simulates the competition between two species of bacteria. The first
% bacteria, the controller, produces antibiotic to slow the growth of the
% second, the invasive species. The problem is non-trivial because
% producing antibiotic puts a metabolic load on the controller, slowing its
% own growth rate. The system is modelled by the following system of ODEs
%   R'   = F(R, A)
%   P'wt = G(R, A)
%   A'   = H(A, \mu)
% where R is the ratio of the size of the invasive species population to
% the size of the controller population, Pwt is the ratio of the size of
% the invasive species population comparted to its starting value, and A is
% concentration of the antibiotic. \mu takes on values in the set {0,1}
% where \mu = 0 corresponds to antibiotic not being produced and \mu = 1
% corresponds to the antibiotic being produced.
% 
%The inputs for system are:
% R_0 = the intial ratio of the two population sizes Pwt_0/Pc_0
% rho = the metabolic cost on the controller for producing the antibiotic.
%       values should be in the set [0,100], where 0 corresponds to no cost
%       and larger values correspond to higher costs.
%
%The output of the simulation is a plot of the paramters R, Pwt, A, \mu,
% and the reference size as they vary through time.
%% Initial Conditions and Important Parameters

% Here we choose the inital population ratio. We set the invasize species
% population ratio (relative to its original size) to 1 if the population
% ratio is not zero. If the population ratio is zero, then the invasive
% species population ratio is set to zero, and the simulation is over.
R_ = 1.078749;
if R_ > 0
    PropWT_ = 1;
else
    PropWT_ = 0;
end

% Here the antibiotic concentration is set to zero because before the
% simulation starts, the controller is not producing antibiotic. x is a
% vector containing the starting state of the system.
A_ = 0;
x = [R_ PropWT_ A_ 1];

% The simulatio starts at time zero. The reference size is set to the
% starting size of the invasive species population. Solutions is a vector
% corresponding to the parameters that will be plotted. Note that solutions
% will grow to n x 5 matrix as the simulation runs.
start_time = 0;
WT_ref     = PropWT_;
solutions  = [start_time R_ PropWT_ A_ 1];

% recession_length corresponds to how far time steps ahead on each
% iteration. Higher step counts increase resolution. newIndex corresponds
% to time step that is recessionLength minutes into the simulation.
lookahead = 30;
recessionLength = 10;
ratio = recessionLength/lookahead;
stepCount = 100;
newIndex = cast(stepCount*ratio, 'int32');
rho = .01;
endSimulation = 10000;
%% Perform RHC

iteration = 0;

while iteration < endSimulation/recessionLength
    
    % Generates the each timestep in the iteration.
    end_time = start_time + lookahead;
    tspan=start_time:lookahead/stepCount:end_time;

    solution = DifferentialEquations(x,tspan,rho);
    
    x = [solution(newIndex+1,1), ...
         solution(newIndex+1,2), ...
         solution(newIndex+1,3), ...
         1];
     
    solutions = [solutions; ...
                 tspan(2:newIndex+1).', ...
                 solution(2:newIndex+1,:)];
    
    start_time = start_time + recessionLength;
    
    iteration = iteration + 1;
end
%% Plot Solutions

% Plots the population ratio in red.
plot(solutions(:,1),(solutions(:,2)),'LineWidth',2,'Color',[1 0 0])
hold on
% Plots the invasive population ratio in blue.
plot(solutions(:,1),(solutions(:,3)),'LineWidth',2,'Color',[0 0 1])
hold on
% Plots the antibiotic concentraion in a green.
plot(solutions(:,1),.01*(solutions(:,4)),'LineWidth',2,'Color',[0 1 0])
hold on
% Plots the controller state in a seafoam.
plot(solutions(:,1),0.1*R_*(solutions(:,5)),'LineWidth',2,'Color',[0.5 1 0])

legend('Population Ratios (Wild-Type/Controller)', ...
       'Wild-Type Proportion (Current/Initial)', ...
       'Antiobiotic Concentration', ...
       '\mu')
   
hold off

xlim([0 0.1*endSimulation])
xlabel('Time (minutes)')
title('Population Dynamics')
grid on
