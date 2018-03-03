clear 
close all
% clc

% Problem settings
Parameters = struct();
Parameters.DataSet = 'Heart.mat';
% Parameters.DataSet = 'PimaIndiansDiabetes';
% Parameters.DataSet = 'AussieCredit.mat';
% Parameters.DataSet = 'GermanCredit.mat';
load(Parameters.DataSet,'Like')
Parameters.Dim = size(Like.ModelMatrix,2);

% Algorithmic setting
N = 2^10;
M = 1;
T = 20;
Stepsize = 0.01;
Parameters.Particles = N;
Parameters.Steps = T;
Parameters.TerminalTime = Stepsize * T;
Parameters.MCMCmoves = M;
    
TotalTime = tic;
    SMC = AIS_Resample(Parameters); 
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
    ylabel('$\log\,G_t$','FontSize',25,'Interpreter','LaTeX')
