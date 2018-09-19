%% Initial Conditions 
WT_0 = 200;
C_0 = 50;
A_0 = 0;
x1 = [WT_0, C_0, A_0, 0];
x2 = [WT_0, C_0, A_0, 1];

start_time = 0;
WT_ref     = WT_0;
solutions  = [start_time, WT_0, C_0, A_0, 0, WT_ref];
distances  = cast(intmax,'double');

lookahead  = 30;
recession_length = 10;
ratio = recession_length/lookahead;
step_count = 100;
new_index = cast(step_count*ratio, 'int32');
rho        = .004;
bound      = 250;

%% Perform RHC

count = 0;

while count < 1000
    
    end_time = start_time + lookahead;
    tspan=start_time:lookahead/step_count:end_time;
    WT_ref_vec = WT_ref.*ones(size(tspan)).';
    
    [times1, solutions1] = simulate(x1,tspan,rho);
    [times2, solutions2] = simulate(x2,tspan,rho);
    distance1 = calculate_distance(WT_ref_vec,solutions1(:,1));
    distance2 = calculate_distance(WT_ref_vec,solutions2(:,1));
    
    if distance1 <= distance2
        x1 = [solutions1(new_index+1,1),solutions1(new_index+1,2), ...
              solutions1(new_index+1,3),0];
        x2 = [solutions1(new_index+1,1),solutions1(new_index+1,2), ...
              solutions1(new_index+1,3),1];
        solutions = [solutions; tspan(2:new_index+1).', ...
                     solutions1(2:new_index+1,:), WT_ref_vec(2:new_index+1)];
        distances = [distances;distance1];
        start_time = start_time + 10;
    else
        x1 = [solutions2(new_index+1,1),solutions2(new_index+1,2), ...
              solutions2(new_index+1,3),0];
        x2 = [solutions2(new_index+1,1),solutions2(new_index+1,2), ...
              solutions2(new_index+1,3),1];
        solutions = [solutions; tspan(2:new_index+1).', ...
                     solutions2(2:new_index+1,:), WT_ref_vec(2:new_index+1)];
        distances = [distances;distance2];
        start_time = start_time + recession_length;
    end
    
    
    if distances(end) < bound && distances(end-1) < bound
        WT_ref = .99*WT_ref;
    end
    
    
    count = count + 1;
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
axis([0 1000 0 400])
xlabel('Time (minutes)')
title('Population Dynamics')
grid on

