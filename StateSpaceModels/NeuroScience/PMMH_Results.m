clear 
close all
clc

% Autocorrelation of PMMH chains 
load('PMMH_BPF.mat')
load('PMMH_cSMC.mat')

nlags = 50;
lags = 0:nlags;
        
% Compute BPF autocorrelation
[BPF_alpha_acf,~,~] = autocorr(theta_BPF(:,1),nlags);
[BPF_sigmasq_acf,~,~] = autocorr(theta_BPF(:,2),nlags);

% Compute cSMC autocorrelation
[cSMC_alpha_acf,~,~] = autocorr(theta_cSMC(:,1),nlags);
[cSMC_sigmasq_acf,~,~] = autocorr(theta_cSMC(:,2),nlags);

% Compare autocorrelation
figure
    hold on
    plot(lags, BPF_alpha_acf, 'bx')
    for ilag = 0:nlags
        plot(ilag, cSMC_alpha_acf(ilag+1), 'r.','MarkerSize',15)
        plot([ilag ilag], [0 cSMC_alpha_acf(ilag+1)], 'r-')
    end
    set(gca,'FontSize',15) 
    title('Parameter $\alpha$','FontSize',25,'Interpreter','LaTeX')
    xlabel('$\mathrm{Lags}$','FontSize',25,'Interpreter','LaTeX') 
    ylabel('$\mathrm{ACF}$','FontSize',25,'Interpreter','LaTeX') 
    h = legend('BPF','cSMC','Location','Best');
    set(h,'FontSize',16)
    axis('tight')

figure
    hold on
    plot(lags, BPF_sigmasq_acf, 'bx')
    for ilag = 0:nlags
        plot(ilag, cSMC_sigmasq_acf(ilag+1), 'r.','MarkerSize',15)
        plot([ilag ilag], [0 cSMC_sigmasq_acf(ilag+1)], 'r-')
    end
    set(gca,'FontSize',15) 
    title('Parameter $\sigma^2$','FontSize',25,'Interpreter','LaTeX')
    xlabel('$\mathrm{Lags}$','FontSize',25,'Interpreter','LaTeX') 
    ylabel('$\mathrm{ACF}$','FontSize',25,'Interpreter','LaTeX') 
    h = legend('BPF','cSMC','Location','Best');
    set(h,'FontSize',16)
    axis('tight')
    
%% Effective sample size of PMMH chains
load('PMMH_BPF.mat')
load('PMMH_cSMC.mat')

M = 100000;
nlags = 1000;
lags = 0:nlags;
        
% Compute BPF autocorrelation
[BPF_alpha_acf,~,~] = autocorr(theta_BPF(:,1),nlags);
[BPF_sigmasq_acf,~,~] = autocorr(theta_BPF(:,2),nlags);

% Compute cSMC autocorrelation
[cSMC_alpha_acf,~,~] = autocorr(theta_cSMC(:,1),nlags);
[cSMC_sigmasq_acf,~,~] = autocorr(theta_cSMC(:,2),nlags);

% ESS
BPF_alpha_ICAT = 1 + 2 * sum(BPF_alpha_acf(2:end));
cSMC_alpha_IACT = 1 + 2 * sum(cSMC_alpha_acf(2:end));

BPF_sigmasq_IACT = 1 + 2 * sum(BPF_sigmasq_acf(2:end));
cSMC_sigmasq_IACT = 1 + 2 * sum(cSMC_sigmasq_acf(2:end));

disp(['ESS (and %) for PMMH-BPF alpha chain ' num2str(M / BPF_alpha_ICAT) ' (' num2str(100 / BPF_alpha_ICAT) '%)'])
disp(['ESS (and %) for PMMH-cSMC alpha chain ' num2str(M / cSMC_alpha_IACT) ' (' num2str(100 / cSMC_alpha_IACT) '%)'])    

disp(['ESS (and %) for PMMH-BPF sigmasq chain ' num2str(M / BPF_sigmasq_IACT) ' (' num2str(100 / BPF_sigmasq_IACT) '%)'])
disp(['ESS (and %) for PMMH-cSMC sigmasq chain ' num2str(M / cSMC_sigmasq_IACT) ' (' num2str(100 / cSMC_sigmasq_IACT) '%)'])    

%% Kernel density estimates 
figure
    hold on
    [BPF_alpha_pdf,BPF_alpha_x] = ksdensity(theta_BPF(:,1),'support',[0 1]);
    [cSMC_alpha_pdf,cSMC_alpha_x] = ksdensity(theta_cSMC(:,1),'support',[0 1]);
    plot(BPF_alpha_x,BPF_alpha_pdf,'b-')
    plot(cSMC_alpha_x,cSMC_alpha_pdf,'r-')
    set(gca,'FontSize',15) 
    xlabel('$\alpha$','FontSize',25,'Interpreter','LaTeX')
    ylabel('$\mathrm{Posterior\,\,density}$','FontSize',25,'Interpreter','LaTeX') 
    h = legend('BPF','cSMC');
    set(h,'FontSize',16) 
    hold off

figure
    hold on
    [BPF_sigmasq_pdf,BPF_sigmasq_x] = ksdensity(theta_BPF(:,2));
    [cSMC_sigmasq_pdf,cSMC_sigmasq_x] = ksdensity(theta_cSMC(:,2));
    plot(BPF_sigmasq_x,BPF_sigmasq_pdf,'b-')
    plot(cSMC_sigmasq_x,cSMC_sigmasq_pdf,'r-')
    set(gca,'FontSize',15) 
    xlabel('$\sigma^2$','FontSize',25,'Interpreter','LaTeX')
    ylabel('$\mathrm{Posterior\,\,density}$','FontSize',25,'Interpreter','LaTeX') 
    h = legend('BPF','cSMC');
    set(h,'FontSize',16) 
    hold off
