function y = LogUniformDensity(x)
% Log density of a uniform distribution on [-1,1]
% Input argument: x (1 x 1)  
% Output argument: y (1 x 1)

    if and(x >= -1, x <= 1)
        y = 0;
    else
        y = -realmax;
    end

end

