%--------------------------------------------------------------------------

% initial conditions
WT_0 = 800;
C_0 = 500;
A_0 = 0;
mu = 0;
x0 = [WT_0, C_0, A_0, mu].';

% simulation horizon
TSPAN=[0 100];
JSPAN = [0 20];

% rule for jumps
% rule = 1 -> priority for jumps
% rule = 2 -> priority for flows
rule = 1;

options = odeset('RelTol',1e-6,'MaxStep',.1);

% simulate
[t,j,x] = HyEQsolver(@flow_1,@jump_1,@flow_set_1,@jump_set_1,...
    x0,TSPAN,JSPAN,rule,options,'ode23t');

% plot solution
figure(1) % position
clf
subplot(4,1,1), plotHarc(t,j,x(:,1));
grid on
ylabel('WT')
subplot(4,1,2), plotHarc(t,j,x(:,2));
grid on
ylabel('C')
subplot(4,1,3), plotHarc(t,j,x(:,3));
grid on
ylabel('A')
subplot(4,1,4), plotHarc(t,j,x(:,4));
grid on
ylabel('mu')

