% Examine variance of normalizing constant estimator on neuroscience example
clear
close 
clc

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

% cSMC parameters
N = 2^7;
cSMC_N = N;
I = 3;
cSMC_I = I;
Parameters = struct();
Parameters.Particles = N;
Parameters.Dimension = d;
Parameters.Time = T;
Parameters.Iterations = I;
Parameters.Approximator = 'quadratic';
Purpose = 'loglikelihood';
% Grid sigmasq space with alpha = 1
alpha = 0.99;
nrep = 100;
ngrid = 20;
sigmasqgrid = linspace(0.01,0.2,ngrid);
cSMC_LogNormConst = zeros(nrep,ngrid);

tic
for igrid = 1:ngrid
    igrid
    for irep = 1:nrep
        irep
        cSMC_LogNormConst(irep,igrid) = cSMC_Resample(Parameters, Initial, ConstructTransition([alpha sigmasqgrid(igrid)]), Observation, Purpose); 
    end
end    
toc

save('NormConst_cSMC.mat','cSMC_LogNormConst','cSMC_N','cSMC_I')
    

