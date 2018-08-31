function xplus = jump_1(x)

% values
gamma = log(2)/20;
alpha = 1;

% state
WT = x(1);
C = x(2);
A = x(3);
mu = x(4);

if mu == 0
    xplus = [WT, C, A, 1].';
else
    xplus = [WT, C, A, 0].';
end
end