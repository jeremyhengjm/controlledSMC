% PMMH with BPF on neuroscience example
clear
close 
clc

d = 1; % dimension of state
dtheta = 2; % dimension of parameters theta
load('Data_NeuroScience.mat') % load data
T = length(y) - 1; % length of time series

% MCMC parameters
M = 100000; % MCMC iterations
proposal_std = [0.002 0.01]; % proposal's standard deviation

% Initial distribution 
Initial = struct();
Initial.Mean = 0;
Initial.Cov = 1;
Initial.Precision = 1;
Initial.MeanPrecision = Initial.Mean * Initial.Precision;
Initial.LogDetCov = log(Initial.Cov);
Initial.QuadForm = Initial.MeanPrecision * Initial.Mean;
Initial.Sample = @(n) mvnrnd(Initial.Mean, Initial.Cov, n); 

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

% Prior distribution on theta (theta(1) = alpha, theta(2) = sigmasq)
PriorLogDensity = @(x) LogUniformDensity(x(1)) + LogInverseGamma(x(2), 1, 0.1);

% BPF parameters
N = 5529; % cSMC uses I = 3
I = 0;
Parameters = struct();
Parameters.Particles = N;
Parameters.Dimension = d;
Parameters.Time = T;
Parameters.Iterations = I;
Parameters.Approximator = 'quadratic';
Purpose = 'loglikelihood';

tic
% Initialise MCMC
theta_BPF = zeros(M+1,dtheta);
current_theta = [0.99 0.11]; % near MLE
current_logPriorDensity = PriorLogDensity(current_theta);
theta_BPF(1,:) = current_theta;
current_logNormConst = cSMC_Resample(Parameters, Initial, ConstructTransition(current_theta), Observation, Purpose); 
naccept = 0; 

for m = 1:M
    m
    if (mod(m,1000) == 0)
        save('PMMH_BPF.mat','theta_BPF', 'naccept', 'm')
    end
    proposal_theta = current_theta + proposal_std .* randn(1,dtheta);
    proposal_logPriorDensity = PriorLogDensity(proposal_theta);
    if and(abs(proposal_theta(1)) <= 1 , proposal_theta(2) > 0)
        proposal_logNormConst = cSMC_Resample(Parameters, Initial, ConstructTransition(proposal_theta), Observation, Purpose); 
    else
        proposal_logNormConst = -realmax;
        disp('Outside support')
    end
    LogMetropolisRatio = proposal_logPriorDensity + proposal_logNormConst - current_logPriorDensity - current_logNormConst;
    LogAcceptanceProb = LogMetropolisRatio * (LogMetropolisRatio < 0);
    if ( log(rand(1)) < LogAcceptanceProb )
        naccept = naccept + 1;
        current_theta = proposal_theta;
        current_logPriorDensity = proposal_logPriorDensity;
        current_logNormConst = proposal_logNormConst;
    end
    theta_BPF(m+1,:) = current_theta;
end
toc

save('PMMH_BPF.mat','theta_BPF', 'naccept', 'm')