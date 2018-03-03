function y = LogInverseGamma(x,alpha,beta)
% Log density of Inverse Gamma distribution
% Input argument: x (1 x 1)
%                 alpha (1 x 1)
%                 beta (1 x 1)
% Output argument: y (1 x 1)
    
    if (x > 0)        
        y = - (alpha + 1) * log(x) - beta ./ x;
    else
        y = -realmax;
    end


end

