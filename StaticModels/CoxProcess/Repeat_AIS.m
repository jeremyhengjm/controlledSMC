clear
close 
clc
    
% Problem specification
Grid = 30;
d = Grid^2;
    
% Algorithmic settings
Parameters = struct();
N = 2^12 * 5; M = 1; T = 20; StepSize = 0.5;
AIS_N = N;

Parameters.Dim = d;
Parameters.Grid = Grid;
Parameters.Particles = N;
Parameters.MCMCmoves = M;
Parameters.Steps = T;
Parameters.TerminalTime = T * StepSize;

% Data
load('FinPine.mat')

% Call AIS
nreps = 100;
AIS_ESS = zeros(nreps,T+1);
AIS_LogNormConst = zeros(nreps,1);
tic
for irep = 1:nreps
    irep
    AIS = AIS_Resample(Parameters,Like); 
    AIS_ESS(irep,:) = AIS.ESS;
    AIS_LogNormConst(irep) = AIS.LogNormConst(end);
end
toc

save('Results_AIS.mat','AIS_N','AIS_ESS','AIS_LogNormConst')