%MAIN Simulates the competition between two species of bacteria. The first
% bacteria, the controller, produces antibiotic to slow the growth of the
% second, the invasive species. The problem is non-trivial because
% producing antibiotic puts a metabolic load on the controller, slowing its
% own growth rate. The system is modelled by the following system of ODEs:
%   P'wt = F(Pwt, Pc, A)
%   P'c  = G(Pc, A)
%   A'   = H(A, \mu)
% where Pwt is the population of the invasive species, Pc is the population
% of the controller species, and A is concentration of antibiotic. \mu
% takes on values in the set {0,1} where \mu = 0 corresponds to antibiotic
% not being produced and \mu = 1 corresponds to the antibiotic being
% produced.

%The inputs for system are:
% wT_0 = the intial size of the invasive population.
% C_0  = the inital size of the controller population.
% rho  = the metabolic cost on the controller for producing the antibiotic.
%        values should be in the set [0,100], where 0 corresponds to no
%        cost and larger values correspond to higher costs.
%% Inputs section

WT_0 = 15;
C_0  = 5;
rho  = .01;
%% Other variables/constants

A_0 = 0;
x1 = [WT_0, C_0, A_0, 1];

start_time = 0;
WT_ref     = WT_0;
solutions  = [start_time, WT_0, C_0, A_0, 0, WT_ref];
distances  = cast(intmax,'double');

lookahead  = 30;
recessionLength = 10;
ratio = recessionLength/lookahead;
step_count = 100;
new_index = cast(step_count*ratio, 'int32');
threshold      = 10000000000000;
endSimulation  = 10000;
%% Perform RHC

iteration = 0;

while iteration < endSimulation/recessionLength
    
    % Generates the each timestep in the iteration.
    end_time = start_time + lookahead;
    tspan=start_time:lookahead/step_count:end_time;
    
    % Sets the reference for each time in the time span.
    WT_ref_vec = WT_ref.*ones(size(tspan)).';
    
    % Simulates for \mu = 0 and \mu = 1. The size of the invasive species
    % is then compared to the reference size at each time step, and a
    % distance is calculated for both values of \mu.
    [times1, solutions1] = differential_equations(x1,tspan,rho);
    distance1 = calculate_distance(WT_ref_vec,solutions1(:,1));
    
    % If \mu = 0 gets the size of the invasive population closer to the
    % reference, then the controller chooses to have the antibiotic off for
    % the next recessionLength minutes of the simulation. Solutions gets
    % newIndex new "snapshots" of the parameters appended to it. Time moves
    % forward recessionLength minutes and the initial conditions are reset.
        x1 = [solutions1(new_index+1,1), ...
              solutions1(new_index+1,2), ...
              solutions1(new_index+1,3), ...
              1];

        solutions = [solutions; ...
                     tspan(2:new_index+1).', ...
                     solutions1(2:new_index+1,:), ...
                     WT_ref_vec(2:new_index+1)];
        distances = [distances;distance1];
        start_time = start_time + recessionLength;
        
    % If two succesive iterations are small, then the reference size is
    % reduced to 95 percent of its value.
    if distances(end) < threshold && distances(end-1) < threshold
        WT_ref = .95*WT_ref;
    end
    
    
    iteration = iteration + 1;
end
%% Plot Solutions

plot(solutions(:,1),(solutions(:,2)),'LineWidth',2,'Color',[1 0 0])
hold on
plot(solutions(:,1),(solutions(:,3)),'LineWidth',2,'Color',[0 0 1])
hold on
plot(solutions(:,1),(solutions(:,6)),'LineWidth',2,'Color',[1 0.5 0])
hold on
plot(solutions(:,1),0.5*C_0*(solutions(:,5)),'LineWidth',1.5,'Color',[0.5 1 0])
hold on
plot(solutions(:,1),(solutions(:,4)),'LineWidth',2,'Color',[0 1 0])
legend('Invasive Pop.', 'Controller Pop.', 'Target Size', 'State','A')
hold off
xlim([0 100])
xlabel('Time (minutes)')
title('Population Dynamics')
grid on

