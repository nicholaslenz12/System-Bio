%% PLOT_SIMUATION

addpath('../Models')

%% ------------------------------------------------------------------------
% Inputs
% -------------------------------------------------------------------------
WT_0 = 100;
C_0  = 50;
A_0 = 0;

%% ------------------------------------------------------------------------
% SIMULATION PARAMETERS
% -------------------------------------------------------------------------
lookahead = 20;
recessionLength = 10;
step_count = 200;
endSimulation  = 200;
ratio = recessionLength/lookahead;
new_index = cast(step_count*ratio, 'int32');
threshold = 20;
distances = threshold + .01;

%% ------------------------------------------------------------------------
% INITIALIZATION
% -------------------------------------------------------------------------
start_time = 0;
x1     = [WT_0, C_0, A_0, 1];
x2     = [WT_0, C_0, A_0, 2];
x3     = [WT_0, C_0, A_0, 3];
x      = [x1; x2; x3];
scenario_count = size(x(:,1));
WT_ref = 2*WT_0;
solutions = [start_time, WT_0, C_0, A_0, 1, WT_ref];

%% ------------------------------------------------------------------------
% LOOP
% -------------------------------------------------------------------------
iteration = 0;

while iteration < endSimulation/recessionLength

    % Generates the each timestep in the iteration.
    end_time = start_time + lookahead;
    tspan=start_time:lookahead/step_count:end_time;

    % Sets the reference for each time in the time span.
    WT_ref_vec = WT_ref.*ones(size(tspan)).';

    predicted_solutions = mpc(x, tspan);
    distances = zeros(scenario_count);
    for idx=1:scenario_count
        distances(idx) = (sum((predicted_solutions(:,1,idx)-WT_ref_vec).^2))^(1/2);
    end

    [min_distance, min_idx] = min(distances);

    for idx=1:scenario_count
        x(idx,:) = [predicted_solutions(new_index+1,1:end-1,min_idx) idx];
    end
    
    solutions = [solutions; ...
                tspan(2:new_index+1).', ...
                predicted_solutions(2:new_index+1,:,min_idx), ...
                WT_ref_vec(2:new_index+1)];
%     solutions
%     input('>>>')
    distances = [distances; min_distance];

    start_time = start_time + recessionLength;
    iteration = iteration + 1;
end

%% ------------------------------------------------------------------------
% PLOT SOLUTIONS
% -------------------------------------------------------------------------
master_xlim = [0 endSimulation];
figure('Renderer', 'painters', 'Position', [720 450 600 400])

%% ------------------------------------------------------------------------
% POPULATION SIZE PLOT
% -------------------------------------------------------------------------
subplot(2,3,[1 2 4 5])
plot(solutions(:,1),(solutions(:,6)),'LineWidth',1.5,'Color',[0 1 0])
hold on
plot(solutions(:,1),(solutions(:,2)),'LineWidth',1.5,'Color',[1 0 0])
hold on
plot(solutions(:,1),(solutions(:,3)),'LineWidth',1.5,'Color',[0 0 1])
hold on
xlim(master_xlim)
ylim([0 inf])
legend('Target Pop.', 'Invasive Pop.', 'Controller Pop.')
grid on
xlabel('Time')
ylabel('Population Size')
title('Bacteria Populations')

%% ------------------------------------------------------------------------
% ANTIBIOTIC CONCENTRATION PLOT
% -------------------------------------------------------------------------
subplot(2,3,3)
plot(solutions(:,1),(solutions(:,4)),'LineWidth',1.5,'Color',[1 0.5 0])
hold on
xlim(master_xlim)
grid on
xlabel('Time')
ylabel('Concentration')
title('Antibiotic Concentration')

%% ------------------------------------------------------------------------
% SWITCH STATE PLOT
% -------------------------------------------------------------------------
subplot(2,3,6)
switch_state = area(solutions(:,1),(solutions(:,5)));
set(switch_state,'facealpha',.5)
xlim(master_xlim)
grid on
xlabel('Time')
ylabel('On / Off')
title('Switch State')

function solutions = mpc(x, tspan)
    solutions = [];
    for idx=1:size(x,1)
        solutions = cat(3,solutions,model_1(x(idx,:), tspan));
    end
end