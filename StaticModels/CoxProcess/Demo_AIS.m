clear
close all
clc

% Load dataset
load('FinPine.mat')

% Problem specification
Grid = 30;
d = Grid^2;
Parameters = struct();
Parameters.Dim = d;
Parameters.Grid = Grid;

%% Annealed importance sampler with Metropolis-adjusted Langevin algorithm moves
N = 2^12 * 5; M = 1; T = 20; StepSize = 0.5; 
Parameters.Particles = N;
Parameters.Steps = T;
Parameters.TerminalTime = T * StepSize;
Parameters.MCMCmoves = M;
    
% Run AIS
TotalTime = tic;
    SMC = AIS_Resample(Parameters,Like); 
toc(TotalTime)

%% Effective sample size 
figure
    hold on
    plot(0:T,SMC.ESS / N * 100,'b-*')
    axis([0 T 0 100])
    set(gca,'FontSize',15) 
    xlabel('$t$','FontSize',25,'Interpreter','LaTeX')
    ylabel('$ESS\%$','FontSize',25,'Interpreter','LaTeX') 
    
%% Acceptance probabilities
figure
    hold on
    errorbar(1:T,mean(SMC.AvgAcceptProb),(SMC.AvgAcceptProb(2,:) - SMC.AvgAcceptProb(1,:)) / 2,'b')
    set(gca,'FontSize',15) 
    xlabel('$t$','FontSize',25,'Interpreter','LaTeX')
    ylabel('Acceptance probabilties','FontSize',15)
    axis('tight')    
    
%% Normalising constant estimation
    figure
        hold on
        plot(0:T,SMC.LogNormConst,'b-*','LineWidth',1)
        set(gca,'FontSize',15) 
        xlabel('$t$','FontSize',25,'Interpreter','LaTeX')
        ylabel('$\log Z_t$','FontSize',25,'Interpreter','LaTeX') 
        axis('tight')
        
%% Variance of log incremental weights
figure
    plot(0:T,SMC.LogIncrementalWeight)
    axis('tight')
    set(gca,'FontSize',15) 
    xlabel('$t$','FontSize',25,'Interpreter','LaTeX')
    ylabel('$\log\,w_t(X_t)$','FontSize',25,'Interpreter','LaTeX')
