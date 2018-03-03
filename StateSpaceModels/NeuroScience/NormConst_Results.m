clear
close all
clc

% Relative variance of log normalizing constant estimates
ngrid = 20;
sigmasqgrid = linspace(0.01,0.2,ngrid);
load('NormConst_BPF.mat')    
load('NormConst_cSMC.mat')

figure
    hold on
    plot(sigmasqgrid, log10(var(BPF_LogNormConst) ./ (mean(BPF_LogNormConst).^2)), 'b-*')
    plot(sigmasqgrid, log10(var(cSMC_LogNormConst) ./ (mean(cSMC_LogNormConst).^2)), 'r-*')
    set(gca,'FontSize',15) 
    xlabel('$\sigma^2$','FontSize',25,'Interpreter','LaTeX') 
    ylabel('$\log_{10}(\mathrm{RVAR})$','FontSize',25,'Interpreter','LaTeX')
    h = legend('BPF','cSMC');
    set(h,'FontSize',16)
    axis tight
%     savefig('NormConst')
%     print('NormConst','-depsc')


    