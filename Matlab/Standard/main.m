%% Initial Conditions 
R_0 = 4; % Ratio of the two populations
Prop_WT_0 = 1; % Proportion of the WT population present at any time compared to the start
A_0 = 0; % Antibiotic concentration
x1 = [R_0 Prop_WT_0 A_0 0];
x2 = [R_0 Prop_WT_0 A_0 1];

start_time = 0;
WT_ref     = Prop_WT_0;
solutions  = [start_time R_0 Prop_WT_0 A_0 0 WT_ref];

lookahead  = 30;
recession_length = 10;
ratio = recession_length/lookahead;
step_count = 200;
new_index = cast(step_count*ratio, 'int32');
rho        = .005;
bound      = 100;
distances  = bound + .01;

%% Perform RHC

count = 0;

while count < 1000
    
    % Generates the timespan for the iteration.
    end_time = start_time + lookahead;
    tspan=start_time:lookahead/step_count:end_time;
    
    % Sets the reference for each time in the time span.
    WT_ref_vec = WT_ref.*ones(size(tspan)).';
    
    % Simulates for \mu = 0 and \mu = 1
    [times1, solutions1] = simulate(x1,tspan,rho);
    [times2, solutions2] = simulate(x2,tspan,rho);
    distance1 = calculate_distance(WT_ref_vec,solutions1(:,2));
    distance2 = calculate_distance(WT_ref_vec,solutions2(:,2));
    
    if distance1 <= distance2
        x1 = [solutions1(new_index+1,1),solutions1(new_index+1,2),solutions1(new_index+1,3),0];
        x2 = [solutions1(new_index+1,1),solutions1(new_index+1,2),solutions1(new_index+1,3),1];
        solutions = [solutions; tspan(2:new_index+1).', ...
                     solutions1(2:new_index+1,:), WT_ref_vec(2:new_index+1)];
        distances = [distances;distance1];
        start_time = start_time + recession_length;
    else
        x1 = [solutions2(new_index+1,1),solutions2(new_index+1,2),solutions2(new_index+1,3),0];
        x2 = [solutions2(new_index+1,1),solutions2(new_index+1,2),solutions2(new_index+1,3),1];
        solutions = [solutions; tspan(2:new_index+1).', ...
                     solutions2(2:new_index+1,:), WT_ref_vec(2:new_index+1)];
        distances = [distances;distance2];
        start_time = start_time + recession_length;
    end
    
    if distances(end) < bound && distances(end-1) < bound
        WT_ref = .95*WT_ref;
    end
    
    count = count + 1;
end

%% Plot Solutions
plot(solutions(:,1),(solutions(:,2)),'LineWidth',2,'Color',[1 0 0])
hold on
plot(solutions(:,1),(solutions(:,3)),'LineWidth',2,'Color',[0 0 1])
hold on
plot(solutions(:,1),.01*(solutions(:,4)),'LineWidth',2,'Color',[0.5 1 0])
hold on
plot(solutions(:,1),(solutions(:,5)),'LineWidth',2,'Color',[0 1 0])
hold on
plot(solutions(:,1),(solutions(:,6)),'LineWidth',2,'Color',[0 1 1])
legend('Population Ratios (Wild-Type/Controller)', ...
       'Wild-Type Proportion (Current/Initial)','Antiobiotic Concentration', ...
       '\mu','Reference Wild-Type Proportion')
hold off

xlim([0 1000])
xlabel('Time (minutes)')
title('Population Dynamics')
grid on

