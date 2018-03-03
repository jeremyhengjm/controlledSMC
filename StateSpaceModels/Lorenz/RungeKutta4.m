function Y = RungeKutta4(X,interval,F)
%Performs an iteration of classical RK4 numerical integration scheme
% Input arguments
%   X (N x d)                Current position
%   stepsize (1 x 1)         Step size to be taken 
%   F (function)             Velocity field
%
%
% Output arguments
%   XOUT: (1 x dim(X))              Next position
    
    nsteps = 10;
    stepsize = interval / nsteps;
    
    for n = 1:nsteps
        k1 = F(X);
        k2 = F(X + 0.5 * stepsize * k1);
        k3 = F(X + 0.5 * stepsize * k2);
        k4 = F(X + stepsize * k3);
        X = X + (stepsize/6) * (k1 + 2*k2 + 2*k3 + k4);
    end
    Y = X;
    
end 
