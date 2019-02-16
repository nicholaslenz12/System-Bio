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
%% Inputs section ---------------------------------------------------------

WT_0 = 100;
C_0  = 10;
A_0 = 0;
%% Other variables/constants ----------------------------------------------

gamma = .1;

%---- Simulation Parameters ----
lookahead = 20;
recessionLength = 10;
ratio = recessionLength/lookahead;
step_count = 200;
new_index = cast(step_count*ratio, 'int32');
threshold      = 20;
distances      = threshold + .01;
endSimulation  = 1000;

%---- Initialization ----
start_time = 0;
WT_ref = WT_0;
x1     = [WT_0, C_0, A_0, 0];
x2     = [WT_0, C_0, A_0, 1];
x3     = [WT_0, C_0, A_0, 2];
solutions  = [start_time, WT_0, C_0, A_0, 0, WT_ref];
%% Perform RHC ------------------------------------------------------------

iteration = 0;

while iteration < endSimulation/recessionLength
    
    % Generates the each timestep in the iteration.
    end_time = start_time + lookahead;
    tspan=start_time:lookahead/step_count:end_time;
    
    % Sets the reference for each time in the time span.
    WT_ref_vec = WT_ref.*ones(size(tspan)).';
    
    % Simulates the system for selected values of \alpha. The error between
    % WT_ref and WT over 'lookahead' minutes is computed for each
    % simulation.
    solutions1 = differential_equations(x1,tspan);
    solutions2 = differential_equations(x2,tspan);
    solutions3 = differential_equations(x3,tspan);
    distance1 = calculate_distance(WT_ref_vec,solutions1(:,1));
    distance2 = calculate_distance(WT_ref_vec,solutions2(:,1));
    distance3 = calculate_distance(WT_ref_vec,solutions3(:,1));
    
    % If \mu = 0 gets the size of the invasive population closer to the
    % reference, then the controller chooses to have the antibiotic off for
    % the next recessionLength minutes of the simulation. Solutions gets
    % newIndex new "snapshots" of the parameters appended to it. Time moves
    % forward recessionLength minutes and the initial conditions are reset.
    if distance1 == min([distance1 distance2 distance3])
        x1 = [solutions1(new_index+1,1), ...
              solutions1(new_index+1,2), ...
              solutions1(new_index+1,3), ...
              0];
        x2 = [solutions1(new_index+1,1), ...
              solutions1(new_index+1,2), ...
              solutions1(new_index+1,3), ...
              1];
        x3 = [solutions1(new_index+1,1), ...
              solutions1(new_index+1,2), ...
              solutions1(new_index+1,3), ...
              2];
        solutions = [solutions; ...
                     tspan(2:new_index+1).', ...
                     solutions1(2:new_index+1,:), ...
                     WT_ref_vec(2:new_index+1)];
        distances = [distances;distance1];
        start_time = start_time + recessionLength;
    elseif distance2 == min([distance1 distance2 distance3])
        x1 = [solutions2(new_index+1,1), ...
              solutions2(new_index+1,2), ...
              solutions2(new_index+1,3), ...
              0];
        x2 = [solutions2(new_index+1,1), ...
              solutions2(new_index+1,2), ...
              solutions2(new_index+1,3), ...
              1];
        x3 = [solutions2(new_index+1,1), ...
              solutions2(new_index+1,2), ...
              solutions2(new_index+1,3), ...
              2];
        solutions = [solutions; ...
                     tspan(2:new_index+1).', ...
                     solutions2(2:new_index+1,:), ...
                     WT_ref_vec(2:new_index+1)];
        distances = [distances;distance2];
        start_time = start_time + recessionLength;
    else
        x1 = [solutions3(new_index+1,1), ...
              solutions3(new_index+1,2), ...
              solutions3(new_index+1,3), ...
              0];
        x2 = [solutions3(new_index+1,1), ...
              solutions3(new_index+1,2), ...
              solutions3(new_index+1,3), ...
              1];
        x3 = [solutions3(new_index+1,1), ...
              solutions3(new_index+1,2), ...
              solutions3(new_index+1,3), ...
              2];
        solutions = [solutions; ...
                     tspan(2:new_index+1).', ...
                     solutions3(2:new_index+1,:), ...
                     WT_ref_vec(2:new_index+1)];
        distances = [distances;distance3];
        start_time = start_time + recessionLength;
    end
    
    % If two succesive iterations are small, then the reference size is
    % reduced to 95 percent of its value.
    if distances(end) < threshold && distances(end-1) < threshold
        WT_ref = .95*WT_ref;
    end
    
    
    iteration = iteration + 1;
end
%% Plot Solutions ---------------------------------------------------------

master_xlim = [0 endSimulation];
figure('Renderer', 'painters', 'Position', [720 450 600 400])

%---- Population Size Plot ----
subplot(2,3,[1 2 4 5])
plot(solutions(:,1),(solutions(:,2)),'LineWidth',2,'Color',[1 0 0])
hold on
plot(solutions(:,1),(solutions(:,3)),'LineWidth',2,'Color',[0 0 1])
hold on
plot(solutions(:,1),(solutions(:,6)),'LineWidth',2,'Color',[0 1 0])
hold on
xlim(master_xlim)
ylim([0 inf])
legend('Invasive Pop.', 'Controller Pop.','Target Pop.')
grid on
xlabel('Time')
ylabel('Population Size')
title('Bacteria Populations')

%---- Antibiotic Concentration Plot ----
subplot(2,3,3)
plot(solutions(:,1),(solutions(:,4)),'LineWidth',2,'Color',[1 0.5 0])
hold on
plot(master_xlim,[1/gamma 1/gamma],'LineWidth',1.5,'LineStyle','--',...
                                   'Color',[1 0.5 0])
xlim(master_xlim)
grid on
xlabel('Time')
ylabel('Concentration')
title('Antibiotic Concentration')

%---- Switch State Plot ----
subplot(2,3,6)
switch_state = area(solutions(:,1),(solutions(:,5)));
set(switch_state,'facealpha',.5)
xlim(master_xlim)
grid on
xlabel('Time')
ylabel('On / Off')
title('Switch State')


