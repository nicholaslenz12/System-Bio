function [distance] = calculate_distance(vector1,vector2)
%CALCULATEDISTANCE computes the distance between two vectors, vector1 and
%   vector2 given some metric.

distance = (sum((vector1-vector2).^2))^(1/2);
end

