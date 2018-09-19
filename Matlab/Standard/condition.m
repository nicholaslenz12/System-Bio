function [boolean] = condition(values)
%CONDITION returns a boolean value based on the state of the system. To be
%   useful the if-elseif-else block would be more complex.

if values(1) > 1000
    boolean = TRUE;
else
    boolean = FALSE;
end
end

