clear
close 
clc
    
% Problem specification
Grid = 30;
d = Grid^2;

% Algorithmic settings
Parameters = struct();
N = 2^12; I = 3; T = 20; StepSize = 1e-2*5; 
cSMC_N = N;

Parameters.Dim = d;
Parameters.Grid = Grid;
Parameters.Particles = N;
Parameters.Steps = T;
Parameters.TerminalTime = T * StepSize;
Parameters.Iterations = I; 
Parameters.Approximator = 'purequadratic';

% Load dataset
load('FinPine.mat')

% Run controlled SMC
nreps = 100;
cSMC_ESS = zeros(nreps,T+1);
cSMC_LogNormConst = zeros(nreps,1);

tic
for irep = 1:nreps
    irep
    cSMC = cSMC_Resample(Parameters,Like); 
    cSMC_ESS(irep,:) = cSMC{end,2}.ESS; 
    cSMC_LogNormConst(irep) = cSMC{end,2}.LogNormConst(end);
end
toc

save('Results_cSMC.mat','cSMC_N','cSMC_ESS','cSMC_LogNormConst')
