function xdot = flow_1(x)

% values
rmax = log(2)/20;
B1 = 1.7;
K1 = 10;
rho = 0.01;
gamma = log(2)/20;
gammawt = .1;
alpha = 1;

% state
WT = x(1);
C = x(2);
A = x(3);
mu = x(4);
Aeff = C*A/WT;

% differential equations
xdot = zeros(4,1);

xdot(1) = rmax*(1 - (Aeff^B1)/(K1^B1 + Aeff^B1))*x(1) - gammawt*Aeff*x(1)/(K1 + Aeff);
if -rho*x(3) + rmax > 0
    xdot(2) = -rho*x(3) + rmax;
end
if mu == 1
    xdot(3) = alpha - gamma*A;
else
    xdot(3) = -gamma*A;
end

end