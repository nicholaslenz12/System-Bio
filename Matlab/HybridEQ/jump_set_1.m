function inside = jump_set_1(x)

% values
% gamma = log(2)/20;
% alpha = 1;

% state
WT = x(1);
C = x(2);
A = x(3);
mu = x(4);

if WT > 1000 && mu == 0
    inside = 1;
elseif WT < 900 && mu == 1
    inside = 1;
else
    inside = 0;
end
end