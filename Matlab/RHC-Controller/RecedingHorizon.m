function [bool] = RecedingHorizon(populationRatio, rho)
%RECEDINGHORIZON RecedingHorizon simulates two competing bacteria
% populations. To read the full description go to main.m in this directory.
%The inputs for this function are:
% populationRatio : The intial ratio of the two population sizes
%                   Pwt_0/Pc_0.
% rho             : The metablic cost on the controller for producing
%                   antibiotic.
%The outputs are:
% bool            : Has value 1 if the size of the invasive species
%                   population converges, 0 if not.
%% Initial Conditions and Important Parameters

% Here we set the inital population ratio, 'R'. 'PropWT_' takes value 1
% if the initial wild-type population size is non-zero (i.e 'R' is not
% zero). If the wild-type population size is zero ('R' is zero) then
% 'PropWT_' is set to 0.
R_ = populationRatio;
if R_ > 0
    PropWT_ = 1;
else
    PropWT_ = 0;
end

% Here the antibiotic concentration, 'A', is set to 0 because before the
% simulation starts, the controller is not producing antibiotic. 'x1' and
% 'x2' are both vectors that contain the intial conditions of the ODEs.
% They differ in that x1(end) = 0, corresponding to the antibiotic being
% off at first. Similarly, x2(end) = 1, implying the antiobiotic is on at
% the start.
A_ = 0;
x1 = [R_ PropWT_ A_ 0];
x2 = [R_ PropWT_ A_ 1];

% The simulation start time, 'start_time', begins at 0. The reference size,
% 'WT_ref' is set to 'PropWT_'. 'solutions' is a vector corresponding to
% the parameters that will be plotted. Note that 'solutions' will grow to
% n x 6 matrix as the simulation runs.
start_time = 0;
WT_ref     = PropWT_;
solutions  = [start_time R_ PropWT_ A_ 0 WT_ref];

% 'lookahead' sets how far (in time) the simulation runs for on each
% iteration. 'recession_length' corresponds to how far time steps between
% iterations. Higher values of 'stepCount' increases resolution. 'newIndex'
% corresponds to time step that is 'recessionLength' minutes into the
% simulation. Higher values of 'threshold' lets the invasive species
% population vary more from the reference size.
lookahead  = 30;
recessionLength = 10;
ratio = recessionLength/lookahead;
stepCount = 100;
newIndex = cast(stepCount*ratio, 'int32');
threshold      = 100;
distances  = threshold + .01;
%% Perform RHC

while solutions(end,3) > 10e-2 && solutions(end,3) < 10
    
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





