clear all

% initial condition 
ic=zeros(3,1);
ic(1) = 1000;
ic(2) = 300;
ic(3) = 0;

% simulation end time 
tend=1000;

% simulation time span
tspan=0:tend/1000:tend;

%% run ODE
[time,y]=ode45(@growth_control,tspan,ic);
 
WT=y(:,1);
C=y(:,2);
A=y(:,3);

plot(time,WT,'LineWidth',2)
hold on
plot(time,C,'LineWidth',2)
hold on
plot(time,10*A,'LineWidth',2)
hold off

function dxdt = growth_control(t,x) % function that computes dydt


% define parameters here
rmax = log(2)/20;
B1 = 1.7;
K1 = 10;
rho = 1;
gamma = log(2)/20;
gammawt = .1;
alpha = 1;
%period = 2;
%gate = 0.74*period;


% WT=x(1)
% C=x(2)
% A=x(3)

dxdt = zeros(size(x));

Aeff = x(2)*x(3)/x(1);

% define ODEs here
dxdt(1) = rmax*(1 - (Aeff^B1)/(K1^B1 + Aeff^B1))*x(1) - gammawt*Aeff*x(1) ...
          /(K1 + Aeff);

if -rho*x(3) + rmax > 0
    dxdt(2) = -rho*x(3) + rmax;
else
    dxdt(2) = 0;
end
    
if x(1) > 1000
    dxdt(3) = alpha - gamma*x(3);
else
    dxdt(3) = -gamma*x(3);
end

end

