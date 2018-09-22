h = 500;
xrange = [0 6];
yrange  = [0 7];

%{
x = [];
for i=1:h+1
    value = (i-1)*(xrange(end)-xrange(1))/h;
    append = value*ones(1,h+1);
    x = [x,append];
end

y = [];
for i=1:h+1
    append = yrange(1):(yrange(end)-yrange(1))/h:yrange(end);
    y = [y,append];
end
%}
  
%%
% z = zeros(1,(h+1)^2);
% ^^^^^^^^^^^^^^^^^^^^^
% Don't uncomment this!

%{
for i=1:(h+1)^2
    z(i) = RecedingHorizon(x(i),y(i));
    if mod(i,h) == 0
        disp(sprintf('%-15.2d',i));
    end
end
%}

scatter(x(z == 1), y(z == 1), 8.4, 'b', 'filled', 'square')
hold on;
scatter(x(z == 0), y(z == 0), 8.4, 'r', 'filled', 'square')
xlabel('pop ratio')
ylabel('\rho')