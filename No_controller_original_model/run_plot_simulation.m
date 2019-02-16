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
% WT_0 = the intial size of the invasive population.
% C_0  = the inital size of the controller population.
% A_0  = the inital antibiotic concentration.
%% Inputs section

WT_0 = 0;
C_0  = 25;
A_0 = 0;
x1 = [WT_0, C_0, A_0, 1];
%% Other variables/constants

rmax = log(2)/20;
B1 = 1.7;
K1 = 10;
B2 = 1.7;
K2 = 10;
gamma = log(2)/20;
gammawt = .1;
alpha = 1;

start_time = 0;
solutions  = [start_time, WT_0, C_0, A_0, 0];

lookahead  = 30;
recessionLength = 10;
ratio = recessionLength/lookahead;
step_count = 100;
new_index = cast(step_count*ratio, 'int32');
endSimulation  = 10000;
%% Perform RHC

iterations_run = 0;

while iterations_run < endSimulation/recessionLength
    
    % Generates the each timestep in the iteration.
    end_time = start_time + lookahead;
    tspan=start_time:lookahead/step_count:end_time;
    
    % Simulates for \mu = 0 and \mu = 1. The size of the invasive species
    % is then compared to the reference size at each time step, and a
    % distance is calculated for both values of \mu.
    [times1, solutions1] = differential_equations(x1,tspan);
    
    x1 = [solutions1(new_index+1,1), ...
          solutions1(new_index+1,2), ...
          solutions1(new_index+1,3), ...
          1];

    solutions = [solutions; ...
                 tspan(2:new_index+1).', ...
                 solutions1(2:new_index+1,:)];
    start_time = start_time + recessionLength;
    
    iterations_run = iterations_run + 1;
end
%% Plot Solutions

master_xlim = [0 endSimulation];
figure('Renderer', 'painters', 'Position', [720 450 600 400])

% Population Size Plot
subplot(2,3,[1 2 4 5])
plot(solutions(:,1),(solutions(:,2)),'LineWidth',2,'Color',[1 0 0])
hold on
plot(solutions(:,1),(solutions(:,3)),'LineWidth',2,'Color',[0 0 1])
hold on
xlim(master_xlim)
legend('Invasive Pop.', 'Controller Pop.')
grid on
xlabel('Time')
ylabel('Population Size')
title('Bacteria Populations')

% Antibiotic Concentration Plot
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

% Switch State Plot
subplot(2,3,6)
switch_state = area(solutions(:,1),(solutions(:,5)));
set(switch_state,'facealpha',.5)
xlim(master_xlim)
grid on
xlabel('Time')
ylabel('On / Off')
title('Switch State')


