function [distance] = calculate_norm(vector1,vector2)
distance = (sum((vector1-vector2).^2))^(1/2);
end

