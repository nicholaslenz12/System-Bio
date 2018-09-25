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
%The behavior of the controller is determined by an algorithm embedded in
% the simulation. The goal of the algorithm is to reduce the size of the
% invasive species to near 0. The outline of the algorithm is as follows:
% 1. Choose a reference size for the invasive species. The controller will
%    try to keep the invasive species near this reference size.
% 2. The system is simulated for 30 minutes. The controller, with the
%    knowledge of what the behavior of the invasive species during that
%    time is, chooses to turn on or off antibiotic production for the next
%    10 minutes to try to get the invasive species population close to the
%    reference.
% 3. Time is moved forward 10 minutes. If the invasive species population
%    is close to the reference for the last two iterations, the reference
%    is reduced to 95 percent of the previous size.
%
%The output of the simulation is a plot of the paramters R, Pwt, A, \mu,
% and the reference size as they vary through time.
%% Initial Conditions and Important Parameters

% Here we choose the inital population ratio. We set the invasize species
% population ratio (relative to its original size) to 1 if the population
% ratio is not zero. If the population ratio is zero, then the invasive
% species population ratio is set to zero, and the simulation is over.
R_ = 4;
if R_ > 0
    PropWT_ = 1;
else
    PropWT_ = 0;
end

% Here the antibiotic concentration is set to zero because before the
% simulation starts, the controller is not producing antibiotic. x1 and x2
% are both vectors that contain the intial conditions of the ODEs. They
% differ in that x1(end) = 0, corresponding to the antibiotic being on at
% first. Similarly, x2(end) = 1, implying the antiobiotic is off at the
% start.
A_ = 0;
x1 = [R_ PropWT_ A_ 0];
x2 = [R_ PropWT_ A_ 1];

% The simulatio starts at time zero. The reference size is set to the
% starting size of the invasive species population. Solutions is a vector
% corresponding to the parameters that will be plotted. Note that solutions
% will grow to n x 6 matrix as the simulation runs.
start_time = 0;
WT_ref     = PropWT_;
solutions  = [start_time R_ PropWT_ A_ 0 WT_ref];

% Sets how far ahead the controller can observe the behavior of the
% invasive species. recession_length corresponds to how far time steps
% ahead on each iteration. Higher step counts increase resolution. newIndex
% corresponds to time step that is recessionLength minutes into the
% simulation. Higher thresholds let the invasive species population vary
% more from the reference size.
lookahead  = 30;
recessionLength = 10;
ratio = recessionLength/lookahead;
stepCount = 200;
newIndex = cast(stepCount*ratio, 'int32');
rho        = .01;
threshold      = 100;
distances  = threshold + .01;
endSimulation = 10000;
%% Perform RHC

iteration = 0;

while iteration < endSimulation/recessionLength
    
    % Generates the each timestep in the iteration.
    end_time = start_time + lookahead;
    tspan=start_time:lookahead/stepCount:end_time;
    
    % Sets the reference for each time in the time span.
    WT_ref_vec = WT_ref.*ones(size(tspan)).';
    
    % Simulates for \mu = 0 and \mu = 1. The size of the invasive species
    % is then compared to the reference size at each time step, and a
    % distance is calculated for both values of \mu.
    [times1, solutions1] = Simulate(x1,tspan,rho);
    [times2, solutions2] = Simulate(x2,tspan,rho);
    distance1 = CalculateDistance(WT_ref_vec,solutions1(:,2));
    distance2 = CalculateDistance(WT_ref_vec,solutions2(:,2));
    
    % If \mu = 0 gets the size of the invasive population closer to the
    % reference, then the controller chooses to have the antibiotic off for
    % the next recessionLength minutes of the simulation. Solutions gets
    % newIndex new "snapshots" of the parameters appended to it. Time moves
    % forward recessionLength minutes and the initial conditions are reset.
    if distance1 <= distance2
        x1 = [solutions1(newIndex+1,1), ...
              solutions1(newIndex+1,2), ...
              solutions1(newIndex+1,3), ...
              0];
        x2 = [solutions1(newIndex+1,1), ...
              solutions1(newIndex+1,2), ...
              solutions1(newIndex+1,3), ...
              1];
        solutions = [solutions; ...
                     tspan(2:newIndex+1).', ...
                     solutions1(2:newIndex+1,:), ...
                     WT_ref_vec(2:newIndex+1)];
        distances = [distances;distance1];
        start_time = start_time + recessionLength;
    % If \mu = 1 gets the size of the invasive population closer to the
    % reference, then the controller chooses to have the antibiotic on for
    % the next recessionLength minutes of the simulation. Solutions gets
    % newIndex new "snapshots" of the parameters appended to it. Time moves
    % forward recessionLength minutes and the initial conditions are reset.
    else
        x1 = [solutions2(newIndex+1,1), ...
              solutions2(newIndex+1,2), ...
              solutions2(newIndex+1,3), ...
              0];
        x2 = [solutions2(newIndex+1,1), ...
              solutions2(newIndex+1,2), ...
              solutions2(newIndex+1,3), ...
              1];
        solutions = [solutions; ...
                     tspan(2:newIndex+1).', ...
                     solutions2(2:newIndex+1,:), ...
                     WT_ref_vec(2:newIndex+1)];
        distances = [distances;distance2];
        start_time = start_time + recessionLength;
    end
    
    % If two succesive iterations are small, then the reference size is
    % reduced to 95 percent of its value.
    if distances(end) < threshold && distances(end-1) < threshold
        WT_ref = .95*WT_ref;
    end
    
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
hold on
% Plots the controller state in a orange.
plot(solutions(:,1),(solutions(:,6)),'LineWidth',2,'Color',[1 0.5 0])

legend('Population Ratios (Wild-Type/Controller)', ...
       'Wild-Type Proportion (Current/Initial)', ...
       'Antiobiotic Concentration', ...
       '\mu', ...
       'Reference Wild-Type Proportion')
   
hold off

xlim([0 0.1*endSimulation])
xlabel('Time (minutes)')
title('Population Dynamics')
grid on
