clear
close all
clc

%% Observation noise = 0.1
load('Results_SigmaSq_point1.mat')
cSMC_N = 2^9;
BPF_N = 2252;

disp('Log Normalising constant')
disp('BPF')
mean(BPF_LogNormConst)
std(BPF_LogNormConst)
disp('cSMC')
mean(cSMC_LogNormConst)
std(cSMC_LogNormConst)

disp('Variance BPF')
var(BPF_LogNormConst)
disp('Variance cSMC')
var(cSMC_LogNormConst)
disp('Variance gain')
var(BPF_LogNormConst) /  var(cSMC_LogNormConst)

disp('ESS')
mean(mean(BPF_ESS * 100 / BPF_N))
mean(mean(cSMC_ESS * 100 / cSMC_N))

disp('Relative variance')
RVAR_BPF = var(BPF_LogNormConst) ./  mean(BPF_LogNormConst)^2
RVAR_cSMC = var(cSMC_LogNormConst) ./  mean(cSMC_LogNormConst)^2
RVAR_BPF / RVAR_cSMC

%% Observation noise = 0.5
load('Results_SigmaSq_point5.mat')
cSMC_N = 2^9;
BPF_N = 4710;

disp('Log Normalising constant')
disp('BPF')
mean(BPF_LogNormConst)
std(BPF_LogNormConst)
disp('cSMC')
mean(cSMC_LogNormConst)
std(cSMC_LogNormConst)

disp('Variance BPF')
var(BPF_LogNormConst)
disp('Variance cSMC')
var(cSMC_LogNormConst)
disp('Variance gain')
var(BPF_LogNormConst) /  var(cSMC_LogNormConst)

disp('ESS')
mean(mean(BPF_ESS * 100 / BPF_N))
mean(mean(cSMC_ESS * 100 / cSMC_N))

disp('Relative variance')
RVAR_BPF = var(BPF_LogNormConst) ./  mean(BPF_LogNormConst)^2
RVAR_cSMC = var(cSMC_LogNormConst) ./  mean(cSMC_LogNormConst)^2
RVAR_BPF / RVAR_cSMC

%% Observation noise = 1
load('Results_SigmaSq_1.mat')
cSMC_N = 2^9;
BPF_N = 6553;

disp('Log Normalising constant')
disp('BPF')
mean(BPF_LogNormConst)
std(BPF_LogNormConst)
disp('cSMC')
mean(cSMC_LogNormConst)
std(cSMC_LogNormConst)

disp('Variance BPF')
var(BPF_LogNormConst)
disp('Variance cSMC')
var(cSMC_LogNormConst)
disp('Variance gain')
var(BPF_LogNormConst) /  var(cSMC_LogNormConst)

disp('ESS')
mean(mean(BPF_ESS * 100 / BPF_N))
mean(mean(cSMC_ESS * 100 / cSMC_N))

disp('Relative variance')
RVAR_BPF = var(BPF_LogNormConst) ./  mean(BPF_LogNormConst)^2
RVAR_cSMC = var(cSMC_LogNormConst) ./  mean(cSMC_LogNormConst)^2
RVAR_BPF / RVAR_cSMC

