%BIFURCATIONPLOT plots the region where there are 3 fixed points for the
% model as of January 2019. 
%% Constants

gammaSteps   = 100;
gammaRange   = [0 gammawt+.004];

%% Grid
x = zeros(gammaSteps+1,1);
for i=1:gammaSteps+1
    x(i) = (i-1)*(gammaRange(2)-gammaRange(1))/(gammaSteps);
end

y = zeros(gammaSteps+1,1);

%%
for i=1:(gammaSteps+1)
    
    c = x(i)
    growth_fixed_gammac = @(x) growth_function(x,c);
    y(i) = fzero(growth_fixed_gammac,[10e-10,10])
    
    if mod(i,gammaSteps) == 0
        disp(sprintf('%d',i))
    end
end

plot(x,y,'LineWidth',2,'Color',[1 0 0])
hold on
t = [0 1];
plot(t,[0 0],'--','Color',[0 0 1],'LineWidth',2)
hold off
xlim(gammaRange)
ylim([-.5 9])
ylabel('\it R','Interpreter','tex')
xlabel('\gamma_{C}','Interpreter','tex')
title('R-coordinate of Fixed Points as a Function of \gamma_{C} (Metabolic Load)')
grid on
legend('Saddle Point','(2) Other Fixed Points')

%% Growth Equation

function growth_output = growth_function(R,gammac)
    alpha = 1;
    rmax = log(2)/20;
    K1 = 10;
    B1 = 1.7;
    K2 = 10;
    B2 = 1.7;
    gamma = log(2)/20;
    gammawt = .1;
    growth_output = rmax*(((alpha/gamma)^B2)/(K2^B2 + (alpha/gamma)^B2)) ...
                  - rmax*(((alpha/gamma/R)^B1)/(K1^B1 + (alpha/gamma/R)^B1)) ...
                  - gammawt*(((alpha/gamma/R))/(K1 + (alpha/gamma/R))) ...
                  + gammac;
end

%rmax*(((alpha/gamma/R)^B1)/(K1^B1 + (alpha/gamma/R)^B1))