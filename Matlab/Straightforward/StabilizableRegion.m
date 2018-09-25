%STABILIZABLEREGION plots the stabilizable region in the for choosen values
% of \rho and R. 
%The inputs are:
% rhoSteps     : The number of increments between successive values of \rho.
% rhoRange     : The range of values of \rho to be plotted.
% rSteps       : The number of increments between successive values of R.
% rRange       : The range of values of R to be plotted.
% fileToLoad : The file storing data from a previous run. It is a matrix of
%   size (rhoSteps+1)x(rSteps+1). Values of 1 correspond to stabilizable
%   points, 0 corresponds to unstabilizable points.
% loadPoint  : The index of the point to start computing stabilizable
%   points from. The stabilizable points are stored in an
%   (rhoSteps+1)x(rSteps+1) matrix, where each element has index:
%   (rSteps+1)*row_index + column_index.

rhoSteps   = 1000;
rhoRange   = [0 100];
rSteps     = 500;
rRange     = [0 6];
fileToLoad = 'test';
loadPoint  = 5000;

x = [];
for i=1:rhoSteps+1
    value = (i-1)*(rhoRange(end)-rhoRange(1))/rhoSteps + rhoRange(1);
    append = value*ones(1,rSteps+1);
    x = [x,append];
end

y = [];
for i=1:rhoSteps+1
    append = rRange(1):(rRange(end)-rRange(1))/rSteps:rRange(end);
    y = [y,append];
end
%%


if exist(fileToLoad, 'file') == 2
    z = importdata(fileToLoad);
else
    z = zeros(rhoSteps+1,rSteps+1);
end

%
for i=loadPoint:(rhoSteps+1)*(rSteps+1)
    z(i) = RecedingHorizon(y(i),x(i));
    if mod(i,rhoSteps) == 0
        disp(sprintf('%d',i))
    end
end
%

%%

scatter(x(z == 1), y(z == 1), 8, 'b', 'filled', 'square')
hold on;
scatter(x(z == 0), y(z == 0), 8, 'r', 'filled', 'square')

%xlim(xrange)
%ylim(yrange)
ylabel('R, population ratio')
xlabel('\rho')
hold off