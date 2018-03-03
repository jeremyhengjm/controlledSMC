clear 
close all
clc

% Load nodel and dataset
Model_NeuroScience


%% Run BPF
N = 2^10;
I = 0;
Parameters = struct();
Parameters.Particles = N;
Parameters.Dimension = d;
Parameters.Time = T;    
Parameters.Iterations = I;
Parameters.Approximator = 'quadratic';
Purpose = 'smoothing';

nreps = 100;
distinct_BPF = zeros(1,nreps);

for irep = 1:nreps
    irep
    cSMC = cSMC_Resample(Parameters, Initial, Transition, Observation, Purpose); 
    distinct_BPF(irep) = length(unique(cSMC(:,1)));
end

%% Run cSMC
N = 2^10;
I = 3;
Parameters = struct();
Parameters.Particles = N;
Parameters.Dimension = d;
Parameters.Time = T;    
Parameters.Iterations = I;
Parameters.Approximator = 'quadratic';
Purpose = 'smoothing';

nreps = 100;
distinct_cSMC = zeros(1,nreps);

for irep = 1:nreps
    irep
    cSMC = cSMC_Resample(Parameters, Initial, Transition, Observation, Purpose); 
    distinct_cSMC(irep) = length(unique(cSMC(:,1)));
end

%%
mean(distinct_cSMC) / mean(distinct_BPF) 
