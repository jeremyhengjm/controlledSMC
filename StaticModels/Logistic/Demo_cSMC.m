clear
close all 
clc

% Problem settings
Parameters = struct();    
Parameters.DataSet = 'Heart.mat';
% Parameters.DataSet = 'PimaIndiansDiabetes';
% Parameters.DataSet = 'AussieCredit.mat';
% Parameters.DataSet = 'GermanCredit.mat';
load(Parameters.DataSet,'Like')
Parameters.Dim = size(Like.ModelMatrix,2);
    
% Algorithmic settings
N = 2^10;
I = 5;
T = 20;
StepSize = 1e-4;
Parameters.Particles = N;
Parameters.Iterations = I; 
Parameters.Steps = T;
Parameters.TerminalTime = T * StepSize;
Parameters.Approximator = 'quadratic';
% Parameters.Approximator = 'purequadratic'; % learn only diagonal matrices

TotalTime = tic;
    cSMC = cSMC_Resample(Parameters); 
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
    h = InsertLegend(I);
    set(h,'FontSize',16)
        
%% Variance of log incremental weights
figure
    subplot(1,2,1)
        plot(0:T,cSMC{1,2}.LogIncrementalWeight)
        uncontrolledrange = axis;
        set(gca,'FontSize',15)
        xlabel('$t$','FontSize',25,'Interpreter','LaTeX')
        ylabel('$\log\,G_t$','FontSize',25,'Interpreter','LaTeX')
    subplot(1,2,2)
        plot(0:T,cSMC{end,2}.LogIncrementalWeight)
        axis(uncontrolledrange)
        set(gca,'FontSize',15)
        xlabel('$t$','FontSize',25,'Interpreter','LaTeX')
        ylabel('$\log\,G_t^{\psi^(I)}$','FontSize',25,'Interpreter','LaTeX')
    