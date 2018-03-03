clear
close 
clc
    
% Algorithmic settings
Parameters = struct();
Parameters.DataSet = 'Heart.mat'
N = 2^10; I = 3; T = 20; StepSize = 1e-4;

% Parameters.DataSet = 'PimaIndiansDiabetes'
% N = 2^10; I = 2; T = 20; StepSize = 1e-3;

% Parameters.DataSet = 'AussieCredit.mat'
% N = 2^10; I = 4; T = 20; StepSize = 1e-3;

% Parameters.DataSet = 'GermanCredit.mat'
% N = 2^10; I = 3; T = 20; StepSize = 5e-4;

load(Parameters.DataSet,'Like')
Parameters.Dim = size(Like.ModelMatrix,2);
Parameters.Particles = N;
Parameters.Iterations = I; 
Parameters.Steps = T;
Parameters.TerminalTime = T * StepSize;
Parameters.Approximator = 'quadratic';
% Parameters.Approximator = 'purequadratic';

nreps = 100;
cSMC_ESS = zeros(nreps,T+1);
cSMC_LogNormConst = zeros(nreps,1);
tic
for irep = 1:nreps
    irep
    cSMC_Output = cSMC_Resample(Parameters); 
    cSMC_ESS(irep,:) = cSMC_Output {end,2}.ESS; 
    cSMC_LogNormConst(irep) = cSMC_Output{end,2}.LogNormConst(end);
end
toc

save('Results_cSMC_Heart.mat','cSMC_ESS','cSMC_LogNormConst')
% save('Results_cSMC_PimaIndians.mat','cSMC_ESS','cSMC_LogNormConst')
% save('Results_cSMC_AussieCredit.mat','cSMC_ESS','cSMC_LogNormConst')
% save('Results_cSMC_GermanCredit.mat','cSMC_ESS','cSMC_LogNormConst')