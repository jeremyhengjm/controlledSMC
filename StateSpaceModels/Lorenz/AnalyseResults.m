clear 
close all
clc

load('APF_VaryF.mat')
load('cSMC_VaryF.mat')

% Vary F 
ngrid = 40;
gridF = linspace(2.5,8.5,ngrid);
nrep = 100;

% Plot likelihood
figure
    subplot(1,2,1)
        hold on
        plot([4.8801 4.8801], [min(min(APF_LogNormConst)) max(max(APF_LogNormConst))],'r-')
        plot(gridF,APF_LogNormConst,'x')
        range = axis;
        xlabel('$F$','FontSize',25,'Interpreter','LaTeX') 
        ylabel('$\log\hat{Z}$','FontSize',25,'Interpreter','LaTeX') 
    
    subplot(1,2,2)
        hold on
        plot([4.8801 4.8801], [min(min(APF_LogNormConst)) max(max(APF_LogNormConst))],'r-')
        plot(gridF,cSMC_LogNormConst,'x')
        axis(range)
        xlabel('$F$','FontSize',25,'Interpreter','LaTeX') 
        ylabel('$\log\hat{Z}$','FontSize',25,'Interpreter','LaTeX') 

% Gain in terms of variance
figure
    hold on
    APF_RVAR = var(APF_LogNormConst) ./ mean(APF_LogNormConst).^2;
    cSMC_RVAR = var(cSMC_LogNormConst) ./ mean(cSMC_LogNormConst).^2;
    plot(gridF,log10(APF_RVAR),'b-*')
    plot(gridF,log10(cSMC_RVAR),'r-*')
    set(gca,'FontSize',15) 
    xlabel('$\alpha$','FontSize',25,'Interpreter','LaTeX') 
    ylabel('$\log_{10}(\mathrm{RVAR})$','FontSize',25,'Interpreter','LaTeX')
    h = legend('APF','cSMC');
    set(h,'FontSize',16)
    axis tight
    
% Average iterations taken by cSMC
figure
    plot(gridF,mean(cSMC_Iterations),'b-*')
    set(gca,'FontSize',15) 
    xlabel('$\alpha$','FontSize',25,'Interpreter','LaTeX') 
    ylabel('$\mathrm{Average\,\,iterations}$','FontSize',25,'Interpreter','LaTeX')
    axis tight
    
% Vary observation noise and dimension
nrep = 100;
dim = 2.^(3:6);
ndim = length(dim);
zetasq = 10.^(-4:-2);
nzetasq = length(zetasq);
iplot = 0;
APF_LogNormConst_Array = zeros(ndim,nrep,nzetasq);
cSMC_LogNormConst_Array = zeros(ndim,nrep,nzetasq);

load('VaryZetasq_Dim8.mat')
APF_LogNormConst_Array(1,:,:) = APF_LogNormConst;
cSMC_LogNormConst_Array(1,:,:) = cSMC_LogNormConst;
disp('For dimension = 8')
disp(['Relative variance of APF (log10 scale) is ' num2str(log10(var(APF_LogNormConst) ./ mean(APF_LogNormConst).^2))])
disp(['Relative variance of cSMC (log10 scale) is ' num2str(log10(var(cSMC_LogNormConst) ./ mean(cSMC_LogNormConst).^2))])
disp(['Gain in variance (log10 scale) is ' num2str(log10(var(APF_LogNormConst) ./ var(cSMC_LogNormConst)))])
disp(['RMSE of APF is ' num2str(sqrt(mean(APF_SquaredError)))])
disp(['RMSE of cSMC is ' num2str(sqrt(mean(cSMC_SquaredError)))])

load('VaryZetasq_Dim16.mat')
APF_LogNormConst_Array(2,:,:) = APF_LogNormConst;
cSMC_LogNormConst_Array(2,:,:) = cSMC_LogNormConst;
disp('For dimension = 16')
disp(['Relative variance of APF (log10 scale) is ' num2str(log10(var(APF_LogNormConst) ./ mean(APF_LogNormConst).^2))])
disp(['Relative variance of cSMC (log10 scale) is ' num2str(log10(var(cSMC_LogNormConst) ./ mean(cSMC_LogNormConst).^2))])
disp(['Gain in variance (log10 scale) is ' num2str(log10(var(APF_LogNormConst) ./ var(cSMC_LogNormConst)))])
disp(['RMSE of APF is ' num2str(sqrt(mean(APF_SquaredError)))])
disp(['RMSE of cSMC is ' num2str(sqrt(mean(cSMC_SquaredError)))])

load('VaryZetasq_Dim32.mat')
APF_LogNormConst_Array(3,:,:) = APF_LogNormConst;
cSMC_LogNormConst_Array(3,:,:) = cSMC_LogNormConst;
disp('For dimension = 32')
disp(['Relative variance of APF (log10 scale) is ' num2str(log10(var(APF_LogNormConst) ./ mean(APF_LogNormConst).^2))])
disp(['Relative variance of cSMC (log10 scale) is ' num2str(log10(var(cSMC_LogNormConst) ./ mean(cSMC_LogNormConst).^2))])
disp(['Gain in variance (log10 scale) is ' num2str(log10(var(APF_LogNormConst) ./ var(cSMC_LogNormConst)))])
disp(['RMSE of APF is ' num2str(sqrt(mean(APF_SquaredError)))])
disp(['RMSE of cSMC is ' num2str(sqrt(mean(cSMC_SquaredError)))])

load('VaryZetasq_Dim64.mat')
APF_LogNormConst_Array(4,:,:) = APF_LogNormConst;
cSMC_LogNormConst_Array(4,:,:) = cSMC_LogNormConst;
disp('For dimension = 64')
disp(['Relative variance of APF (log10 scale) is ' num2str(log10(var(APF_LogNormConst) ./ mean(APF_LogNormConst).^2))])
disp(['Relative variance of cSMC (log10 scale) is ' num2str(log10(var(cSMC_LogNormConst) ./ mean(cSMC_LogNormConst).^2))])
disp(['Gain in variance (log10 scale) is ' num2str(log10(var(APF_LogNormConst) ./ var(cSMC_LogNormConst)))])
disp(['RMSE of APF is ' num2str(sqrt(mean(APF_SquaredError)))])
disp(['RMSE of cSMC is ' num2str(sqrt(mean(cSMC_SquaredError)))])

figure
    for idim = 1:ndim
        for izetasq = 1:nzetasq
            iplot = iplot + 1;
            subplot(ndim,nzetasq,iplot)            
            title(['$d = $ ' num2str(dim(idim)) ', $\zeta^2 = $ ' num2str(zetasq(izetasq)) ],'FontSize',15,'Interpreter','LaTeX')
            aboxplot([squeeze(APF_LogNormConst_Array(idim,:,izetasq))' squeeze(cSMC_LogNormConst_Array(idim,:,izetasq))'],'labels',{'APF' 'cSMC'})
        end
    end

    
    