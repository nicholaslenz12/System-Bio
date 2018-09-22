xsteps     = 100;
xrange     = [3.6 4];
ysteps     = 100;
yrange     = [0 0.5];
fileToLoad = '';
loadPoint  = 1;

x = [];
for i=1:xsteps+1
    value = (i-1)*(xrange(end)-xrange(1))/xsteps + xrange(1);
    append = value*ones(1,ysteps+1);
    x = [x,append];
end

y = [];
for i=1:xsteps+1
    append = yrange(1):(yrange(end)-yrange(1))/ysteps:yrange(end);
    y = [y,append];
end
%%

if exist(fileToLoad, 'file') == 2
    z = importdata(fileToLoad);
else
    %z = zeros(xsteps+1,ysteps+1);
end

%{
for i=loadPoint:(xsteps+1)*(ysteps+1)
    z(i) = RecedingHorizon(x(i),y(i));
    if mod(i,xsteps) == 0
        disp(sprintf('%d',i))
    end
end
%}

%%

scatter(x(z == 1), y(z == 1), 30, 'b', 'filled', 'square')
hold on;
scatter(x(z == 0), y(z == 0), 30, 'r', 'filled', 'square')

%xlim(xrange)
%ylim(yrange)
xlabel('R, population ratio')
ylabel('\rho')
hold off