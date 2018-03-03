clear 
close all
clc

load('Results_AIS.mat')
load('Results_cSMC.mat')

load('TrueLogNormConst.mat')

disp('AIS normalizing constant estimation')
disp('Mean')
mean(AIS_LogNormConst)

disp('Standard deviation')
std(AIS_LogNormConst)

disp('Root mean squared error')
AIS_RMSE = sqrt(mean((AIS_LogNormConst - TrueLogNormConst).^2))

disp('cSMC normalizing constant estimation')
disp('Mean')
mean(cSMC_LogNormConst)
disp('Standard deviation')
std(cSMC_LogNormConst)
disp('Root mean squared error')
cSMC_RMSE = sqrt(mean((cSMC_LogNormConst - TrueLogNormConst).^2))

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
    set(findobj(get(h(1), 'parent'), 'type', 'text'), 'fontsize', 15);
    set(gca,'FontSize',15)
    ylabel('$\log Z$','FontSize',25,'Interpreter','LaTeX') 
    