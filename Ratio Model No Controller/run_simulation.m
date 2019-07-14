function [bool] = run_simulation(populationRatio, rho)
%RECEDINGHORIZON RecedingHorizon simulates two competing bacteria
% populations. To read the full description go to Main.m in this directory.
%The inputs for this function are:
% populationRatio : The intial ratio of the two population sizes
%                   Pwt_0/Pc_0.
% rho             : The metablic cost on the controller for producing
%                   antibiotic.
%The outputs are:
% bool            : Has value 1 if the size of the invasive species
%                   population converges, 0 if not.
%% Initial Conditions and Important Parameters

% Here we choose the inital population ratio. We set the invasize species
% population ratio (relative to its original size) to 1 if the population
% ratio is not zero. If the population ratio is zero, then the invasive
% species population ratio is set to zero, and the simulation is over.
R_ = populationRatio;
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
solutions  = [start_time R_ PropWT_ A_ 1];

% recession_length corresponds to how far time steps ahead on each
% iteration. Higher step counts increase resolution. newIndex corresponds
% to time step that is recessionLength minutes into the simulation.
lookahead  = 30;
recessionLength = 10;
ratio = recessionLength/lookahead;
stepCount = 100;
newIndex = cast(stepCount*ratio, 'int32');
%% Perform RHC

while solutions(end,3) > 10e-2 && solutions(end,3) < 10
    
    % Generates the each timestep in the iteration.
    end_time = start_time + lookahead;
    tspan=start_time:lookahead/stepCount:end_time;

    solution = differential_equations(x,tspan,rho);
    
    x = [solution(newIndex+1,1), ...
         solution(newIndex+1,2), ...
         solution(newIndex+1,3), ...
         1];
     
    solutions = [solutions; ...
                 tspan(2:newIndex+1).', ...
                 solution(2:newIndex+1,:)];
    
    start_time = start_time + recessionLength;
end
%% Determine Convergence

% If the simulation ended because the invasive species populatio ratio was
% small.
if solutions(end,3) <= 10e-2
    bool = 1;
% If the simulation ended because the invasive species populatio ratio was
% large.
else
    bool = 0;
end
end
    
    
    
    
    
    