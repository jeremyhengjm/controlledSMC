clear 
close all
clc

load('Results_AIS_Heart.mat')
load('Results_cSMC_Heart.mat')

% load('Results_AIS_PimaIndians.mat')
% load('Results_cSMC_PimaIndians.mat')

% load('Results_AIS_AussieCredit.mat')
% load('Results_cSMC_AussieCredit.mat')

% load('Results_AIS_GermanCredit.mat')
% load('Results_cSMC_GermanCredit.mat')

load('TrueLogNormConst.mat')

switch DataSet
    case 'PimaIndians'
        TrueLogNormConst = TrueLogNormConst(1);
    case 'Heart'
        TrueLogNormConst = TrueLogNormConst(2);
    case 'AussieCredit'
        TrueLogNormConst = TrueLogNormConst(3);
    case 'GermanCredit'
        TrueLogNormConst = TrueLogNormConst(4);
        
end

%% AIS results
disp('AIS normalizing constant estimation')
disp('Mean')
mean(AIS_LogNormConst)
disp('Standard deviation')
std(AIS_LogNormConst)
disp('Variance')
var(AIS_LogNormConst)
disp('Root mean squared error')
AIS_RMSE = sqrt(mean((AIS_LogNormConst - TrueLogNormConst).^2))

%% cSMC results
disp('cSMC normalizing constant estimation')
disp('Mean')
mean(cSMC_LogNormConst)
disp('Standard deviation')
std(cSMC_LogNormConst)
disp('Variance')
var(cSMC_LogNormConst)
disp('Root mean squared error')
cSMC_RMSE = sqrt(mean((cSMC_LogNormConst - TrueLogNormConst).^2))

%% Compare results
disp('Gain in standard deviation over AIS')
std(AIS_LogNormConst) / std(cSMC_LogNormConst) 
disp('Gain in variance over AIS')
var(AIS_LogNormConst) / var(cSMC_LogNormConst) 
disp('Gain in RMSE over AIS')
AIS_RMSE / cSMC_RMSE

%% Average ESS% over time and repetitions
disp('Average ESS% over time and repetitions')
disp('AIS')
mean(mean(AIS_ESS,2)) * 100 / AIS_N 
disp('cSMC')
mean(mean(cSMC_ESS,2)) * 100 / cSMC_N

%% Boxplot of normalizing constant estimates
figure
    algorithm = {'AIS','cSMC'};
    h = boxplot([AIS_LogNormConst cSMC_LogNormConst],algorithm);
    set(findobj(get(h(1), 'parent'), 'type', 'text'), 'fontsize', 13);
    ylabel('Log normalizing constant')
    