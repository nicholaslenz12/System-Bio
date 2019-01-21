function [distance] = calculate_distance(vector1,vector2)
%CALCULATEDISTANCE computes the Eucllidean distance between two vectors,
% vector1 and vector2.
%The inputs for this function are:
% vector1 : The first vector.
% vector2 : The second vector.
%The outputs are:
% distance : The distance between the vectors.

distance = (sum((vector1-vector2).^2))^(1/2);
end

