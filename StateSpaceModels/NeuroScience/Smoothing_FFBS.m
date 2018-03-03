% Forward filtering and backward smoother on neuroscience example
clear
close 
clc

alpha = 0.99; % near MLE
sigmasq = 0.11;
d = 1; % dimension of state
dtheta = 2; % dimension of parameters theta
load('Data_NeuroScience.mat') % load data
T = length(y) - 1; % length of time series

% Initial distribution 
Initial = struct();
Initial.Mean = 0;
Initial.Cov = 1;
Initial.Precision = 1;
Initial.MeanPrecision = Initial.Mean * Initial.Precision;
Initial.LogDetCov = log(Initial.Cov);
Initial.QuadForm = Initial.MeanPrecision * Initial.Mean;
Initial.Sample = @(n) mvnrnd(Initial.Mean, Initial.Cov, n); 

% Transition kernel f
Transition = struct();
Transition.F = alpha;
Transition.Mean = @(t,x) x * Transition.F;
Transition.G = sigmasq;
Transition.Ginv = 1 / sigmasq;
Transition.LogDetG = log(Transition.G);
Transition.Sample = @(x,n) mvnrnd(Transition.Mean(1,x), Transition.G, n);
Transition.Const = -0.5 * log(2 * pi * sigmasq);
Transition.LogDensity = @(x,y) - 0.5 * (y - alpha * x).^2 / sigmasq;

% Observation
Observation = struct();
expit = @(x) exp(x) ./ (1 + exp(x));
ntrials = 50;
Observation.y = y;
Observation.Const = zeros(T+1,1);
for t = 1:(T+1)
    Observation.Const(t) = log( nchoosek(ntrials, Observation.y(t)) );
end
Observation.LogDensity = @(t,x) Observation.Const(t) + Observation.y(t) * log( expit(x) ) ... 
                                + ( ntrials - Observation.y(t) ) * log( 1 - expit(x) );

% SMC parameters
N = 176; 
Parameters = struct();
Parameters.Particles = N;
Parameters.Dimension = d;
Parameters.Time = T;

nreps = 100;
FFBS_SmoothingPaths = zeros(nreps,N,T+1);
FFBS_SmoothingWeights = zeros(nreps,N,T+1);
for irep = 1:nreps
    irep
    cSMC = FFBS(Parameters, Initial, Transition, Observation);
    FFBS_SmoothingPaths(irep,:,:) = cSMC.Trajectory;
    FFBS_SmoothingWeights(irep,:,:) = cSMC.SmoothingWeights;
end

save('Smoothing_FFBS.mat','FFBS_SmoothingPaths','FFBS_SmoothingWeights')
