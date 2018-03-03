clear
close all
clc

% Problem specification
Grid = 30;
d = Grid^2;

% Algorithmic settings
N = 2^12; I = 3; T = 10; StepSize = 1e-2*10; 
Parameters = struct();
Parameters.Dim = d;
Parameters.Grid = Grid;
Parameters.Particles = N;
Parameters.Steps = T;
Parameters.TerminalTime = T * StepSize;
Parameters.Iterations = I; 
Parameters.Approximator = 'purequadratic';

% Load dataset
load('FinPine.mat')

% Run controlled SMC
TotalTime = tic;
    cSMC = cSMC_Resample(Parameters,Like); 
toc(TotalTime)

%% Effective sample size
colors = DefineColors(I);
figure
    hold on
    for p = 1:(I+1)
        plot(0:T,cSMC{p,2}.ESS / N * 100,'-*','Color',colors{p})
    end
    set(gca,'FontSize',15) 
    xlabel('$t$','FontSize',25,'Interpreter','LaTeX')
    ylabel('$ESS\%$','FontSize',25,'Interpreter','LaTeX') 
    axis([0 T 0 100])
    h = InsertLegend(I);
    set(h,'FontSize',16)

%% Normalising constant estimation
figure
    hold on
    for p = 1:(I+1)
        plot(0:T,cSMC{p,2}.LogNormConst,'-*','Color',colors{p})
    end
    set(gca,'FontSize',15) 
    xlabel('$t$','FontSize',25,'Interpreter','LaTeX')
    ylabel('$\log Z_t$','FontSize',25,'Interpreter','LaTeX') 
    axis('tight')
    h = InsertLegend(I);
    set(h,'FontSize',15)

%% Variance of log incremental weights
figure
    subplot(1,2,1)
        plot(0:T,cSMC{1,2}.LogIncrementalWeight)
        uncontrolledrange = axis;
        set(gca,'FontSize',15) 
        xlabel('$t$','FontSize',25,'Interpreter','LaTeX')
        ylabel('$\log\,W_t(X_{0:t})$','FontSize',25,'Interpreter','LaTeX')
    
    subplot(1,2,2)
        plot(0:T,cSMC{end,2}.LogIncrementalWeight)
        axis(uncontrolledrange)
        set(gca,'FontSize',15) 
        xlabel('$t$','FontSize',25,'Interpreter','LaTeX')
        ylabel('$\log\,W_t^{\psi}(X_{0:t})$','FontSize',25,'Interpreter','LaTeX')
    