clear 
close all
clc
%%
% Relative variance of smoothing averages 
T = 2999;
expit = @(x) exp(x) ./ (1 + exp(x));
load('Smoothing_FFBS.mat')
SmoothingMean_FFBS = squeeze(sum(expit(FFBS_SmoothingPaths) .* FFBS_SmoothingWeights, 2));
% load('Smoothing_cSMC.mat')
SmoothingMean_cSMC = squeeze(mean(expit(cSMC_SmoothingPaths), 2));

figure
    hold on
    plot(0:T, log10(var(SmoothingMean_FFBS) ./ (mean(SmoothingMean_FFBS).^2)), 'b-*')
    plot(0:T, log10(var(SmoothingMean_cSMC) ./ (mean(SmoothingMean_cSMC).^2)), 'r-*')
    set(gca,'FontSize',15) % tick size
    xlabel('$t$','FontSize',25,'Interpreter','LaTeX')
    ylabel('$\log_{10}(\mathrm{RVAR})$','FontSize',25,'Interpreter','LaTeX')
    h = legend('FFBS','cSMC');
    set(h,'FontSize',16)
    savefig('Smoothing')
    print('Smoothing','-depsc')

figure
    hold on
    plot(0:T, mean(SmoothingMean_FFBS), 'b-*')
    plot(0:T, mean(SmoothingMean_cSMC), 'r-*')
    set(gca,'FontSize',13) % tick size
    xlabel('$t$','FontSize',25,'Interpreter','LaTeX')
    ylabel('$\mathrm{Smoothing mean}$','FontSize',25,'Interpreter','LaTeX') 
%    
% figure
%     hold on
%     plot(0:T, std(SmoothingMean_FFBS) ./ std(SmoothingMean_cSMC), 'b-*')
%     plot([0 T], [1 1], 'r-')
%     xlabel('$t$','FontSize',25,'Interpreter','LaTeX')
%     ylabel('Ratio of variance','FontSize',25,'Interpreter','LaTeX') 

