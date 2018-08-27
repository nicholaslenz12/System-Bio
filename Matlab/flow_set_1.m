function [value] = flow_set_1(x)

WT = x(1);
C = x(2);
A = x(3);

if ( WT >= 0 && C >= 0 && A >= 0 )
    value = 1;
else
    value = 0;
end
end