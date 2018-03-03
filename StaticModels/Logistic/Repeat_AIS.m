clear
close all
clc
    
% Algorithmic settings
Parameters = struct();
Parameters.DataSet = 'Heart.mat'
N = 1843; M = 1; T = 20; StepSize = 5e-2;

% Parameters.DataSet = 'PimaIndiansDiabetes'
% N = 2^10; M = 1; T = 20; StepSize = 1e-2;

% Parameters.DataSet = 'AussieCredit.mat'
% N = 1843; M = 1; T = 20; StepSize = 3e-2;

% Parameters.DataSet = 'GermanCredit.mat'
% N = 2048; M = 1; T = 20; StepSize = 1e-2;

load(Parameters.DataSet,'Like')
Parameters.Dim = size(Like.ModelMatrix,2);
Parameters.Particles = N;
Parameters.MCMCmoves = M;
Parameters.Steps = T;
Parameters.TerminalTime = T * StepSize;

nreps = 100;
AIS_ESS = zeros(nreps,T+1);
AIS_LogNormConst = zeros(nreps,1);
tic
for irep = 1:nreps
    irep
    AIS_Output = AIS_Resample(Parameters); 
    AIS_ESS(irep,:) = AIS_Output.ESS;
    AIS_LogNormConst(irep) = AIS_Output.LogNormConst(end);
end
toc

save('Results_AIS_Heart.mat','AIS_ESS','AIS_LogNormConst')
% save('Results_AIS_PimaIndians.mat','AIS_ESS','AIS_LogNormConst')
% save('Results_AIS_AussieCredit.mat','AIS_ESS','AIS_LogNormConst')
% save('Results_AIS_GermanCredit.mat','AIS_ESS','AIS_LogNormConst')