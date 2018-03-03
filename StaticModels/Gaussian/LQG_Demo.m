clear
close 
clc

% % Problem settings
Parameters = struct();
Parameters.Dim = 2^5;
Parameters.Rho = 0.8;
Parameters.Obs = 20 * ones(1,Parameters.Dim);
    
% Algorithmic setting
N = 2^11;
T = 10;
StepSize = 0.1;
Parameters.Particles = N;
Parameters.Steps = T;
Parameters.TerminalTime = T * StepSize;
Parameters.Approximator = 'quadratic';
% Parameters.Approximator = 'purequadratic';

tic, cSMC = LQG(Parameters); toc
Exact = cSMC{2,1}; 
        
%% Effective sample size
colors = DefineColors(2);
figure
    hold on
    for p = 1:2
        plot(0:T,cSMC{p,2}.ESS * 100 / N,'-*','MarkerSize',10,'Color',colors{p})
    end
    set(gca,'FontSize',15)
    xlabel('$t$','FontSize',25,'Interpreter','LaTeX')
    ylabel('$ESS\%$','FontSize',25,'Interpreter','LaTeX') 
    axis([0 T 0 100])
    h = legend('Uncontrolled SMC','Optimally controlled SMC');
    set(h,'FontSize',16)
    
%% Normalising constant estimation
figure
    hold on
    plot([0 T],repmat(cSMC{2,2}.V,1,2),'-','LineWidth',2,'Color',colour_lightblue)
    plot(0:T,-Exact.LogNormConst,'r-','LineWidth',2)
    for p = 1:2
        plot(0:T,-cSMC{p,2}.LogNormConst,'-*','MarkerSize',10,'Color',colors{p})
    end
    set(gca,'FontSize',15)
    xlabel('$t$','FontSize',25,'Interpreter','LaTeX')
    ylabel('$-\log Z_t$','FontSize',25,'Interpreter','LaTeX') 
    axis('tight')
    h = legend('Optimal value','True','Uncontrolled SMC','Optimally controlled SMC');
    set(h,'FontSize',16)
              
