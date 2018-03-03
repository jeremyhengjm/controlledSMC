clear 
close all
clc

% Load nodel and dataset
Model_NeuroScience

%% Run controlled SMC
N = 2^10;
I = 2;
Parameters = struct();
Parameters.Particles = N;
Parameters.Dimension = d;
Parameters.Time = T;    
Parameters.Iterations = I;
Parameters.Approximator = 'quadratic';
Purpose = 'demo';

tic, cSMC = cSMC_Resample(Parameters, Initial, Transition, Observation, Purpose); toc

%% Effective sample size
colors = DefineColors(I);
figure
    hold on
    for i = 1:(I+1)
        plot(0:T,(cSMC{i,2}.ESS / N) * 100,'-*','Color',colors{i})
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
    IterationNormConstEnd = zeros(1,I+1);
    for i = 1:(I+1)
        plot(0:T,cSMC{i,2}.LogNormConst,'-*','Color',colors{i})
    end
    set(gca,'FontSize',15) 
    xlabel('$t$','FontSize',25,'Interpreter','LaTeX')
    ylabel('$\log\, Z_t$','FontSize',25,'Interpreter','LaTeX') 
    axis('tight')
    h = InsertLegend(I);
    set(h,'FontSize',15)

%% Variance of log incremental weights
figure
    subplot(1,2,1)
        plot(0:T,cSMC{1,2}.LogIncrementalWeight)
        hold on
        uncontrolledrange = axis;
        set(gca,'FontSize',15) 
        xlabel('$t$','FontSize',25,'Interpreter','LaTeX')
        ylabel('$\log\,W_t(X_{t})$','FontSize',25,'Interpreter','LaTeX')
    subplot(1,2,2)
            plot(0:T,cSMC{end,2}.LogIncrementalWeight)
            axis(uncontrolledrange)
            hold on
            xlabel('$t$','FontSize',25,'Interpreter','LaTeX')
            ylabel('$\log\,W_t^{\psi^*}(X_{t})$','FontSize',25,'Interpreter','LaTeX')

    
%% Estimated coefficient a 
figure
    hold on
    plot(0:T,cSMC{2,1}.A(:,1,1),'.','Color',colors{2},'MarkerSize',15)
        for i = 2:(I)
            plot(0:T,cSMC{i+1,1}.A(:,1,1)-cSMC{i,1}.A(:,1,1),'.','Color',colors{i+1},'MarkerSize',15)
        end
        set(gca,'FontSize',15) 
        xlabel('$t$','FontSize',25,'Interpreter','LaTeX')
        ylabel('$a_t^{i}$','FontSize',25,'Interpreter','LaTeX') 
        h = legend('Iteration 1','Iteration 2','Iteration 3','Location','Best');   
        set(h,'FontSize',16)
        hold off
        
%% Estimated coefficient b
figure
    hold on
    plot(0:T,cSMC{2,1}.b(:,1),'.','Color',colors{2},'MarkerSize',15)
    for i = 2:I
        plot(0:T,cSMC{i+1,1}.b(:,1)-cSMC{i,1}.b(:,1),'.','Color',colors{i+1},'MarkerSize',15)
    end
    set(gca,'FontSize',15) 
    xlabel('$t$','FontSize',25,'Interpreter','LaTeX')
    ylabel('$b_t^{i}$','FontSize',25,'Interpreter','LaTeX') 
    h = legend('Iteration 1','Iteration 2','Iteration 3','Location','Best');   
    set(h,'FontSize',16)
    hold off

%% Estimated coefficient c
figure
    hold on
    plot(0:T,cSMC{2,1}.c,'.','Color',colors{2},'MarkerSize',15)
    for i = 2:I
        plot(0:T,cSMC{i+1,1}.c-cSMC{i,1}.c,'.','Color',colors{i+1},'MarkerSize',15)
    end
    set(gca,'FontSize',15) 
    xlabel('$t$','FontSize',25,'Interpreter','LaTeX')
    ylabel('$c_t^{i}$','FontSize',25,'Interpreter','LaTeX') 
    h = legend('Iteration 1','Iteration 2','Iteration 3','Location','Best');   
    set(h,'FontSize',16)
    hold off
