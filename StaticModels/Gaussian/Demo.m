clear
close all
clc

% Problem settings
Parameters = struct();
Parameters.Demo = 'yes'; % controls whether we do LQG
Parameters.Dim = 2^5;
Parameters.Rho = 0.8;
Parameters.Obs = 20 * ones(1,Parameters.Dim);

% Algorithmic settings
N = 2^11;
T = 10;
I = 2;
StepSize = 0.1;
Parameters.Particles = N;
Parameters.Steps = T;
Parameters.TerminalTime = T * StepSize;
Parameters.Iterations = I; % no of iterations
Parameters.Approximator = 'quadratic';
%     Parameters.Approximator = 'purequadratic'; % learn only diagonal matrices

TotalTime = tic;
    [cSMC,Exact] = cSMC_Resample(Parameters); 
toc(TotalTime)

%% Effective sample size
colors = DefineColors(I);

figure
    hold on
    for p = 1:(I+1)
        plot(0:T,cSMC{p,2}.ESS * 100 / N,'-*','Color',colors{p})
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
    plot(0:Parameters.Steps,Exact.LogNormConst,'r-','LineWidth',3)
    for p = 1:(I+1)
        plot(0:T,cSMC{p,2}.LogNormConst,'-*','Color',colors{p})
    end
    set(gca,'FontSize',15) 
    xlabel('$t$','FontSize',25,'Interpreter','LaTeX')
    ylabel('$\log Z_t$','FontSize',25,'Interpreter','LaTeX') 
    switch I
        case 2
            h = legend('True','Uncontrolled','Iteration 1','Iteration 2','Location','Best'); % p = 2            
        case 4
            h = legend('True','Uncontrolled','Iteration 1','Iteration 2','Iteration 3','Iteration 4','Location','Best'); % p = 4    
        case 6
            h = legend('True','Uncontrolled','Iteration 1','Iteration 2','Iteration 3','Iteration 4', ... 
                        'Iteration 5','Iteration 6','Location','Best'); % p = 6
        case 8
            h = legend('True','Uncontrolled','Iteration 1','Iteration 2','Iteration 3','Iteration 4', ... 
                        'Iteration 5','Iteration 6','Iteration 7','Iteration 8','Location','Best'); % p = 6
    end
    set(h,'FontSize',15)
        

%% Variance of log incremental weights
figure
    subplot(1,2,1)
        plot(0:T,cSMC{1,2}.LogIncrementalWeight)
        hold on
        set(gca,'FontSize',15)
        xlabel('$t$','FontSize',25,'Interpreter','LaTeX')
        ylabel('$\log\,G_t$','FontSize',25,'Interpreter','LaTeX')
        uncontrolledrange = axis;
        
    subplot(1,2,2)
        plot(0:T,cSMC{end,2}.LogIncrementalWeight)
        axis(uncontrolledrange)
        hold on
        set(gca,'FontSize',15)
        xlabel('$t$','FontSize',25,'Interpreter','LaTeX')
        ylabel('$\log\,G_t^{\psi^{(I)}}$','FontSize',25,'Interpreter','LaTeX')
        
    
%% Estimated coefficients
colors = {'b','k','g'};
markersizes = {25,25,12};

figure
    hold on
    plot(0:T,Exact.A(:,1,1),'rx','MarkerSize',20)
    for p = 1:3
        plot(0:T,cSMC{p,1}.A(:,1,1),'.','Color',colors{p},'MarkerSize',markersizes{p})
    end
    set(gca,'FontSize',15) 
    xlabel('$t$','FontSize',25,'Interpreter','LaTeX')
    ylabel('$A_t(1,1)$','FontSize',25,'Interpreter','LaTeX') 
    h = legend('True','Uncontrolled','Iteration 1','Iteration 2','Location','Best');
    set(h,'FontSize',16)
    hold off
    
figure
    hold on
    plot(0:T,Exact.A(:,1,2),'rx','MarkerSize',20)
    for p = 1:3
        plot(0:T,cSMC{p,1}.A(:,1,2),'.','Color',colors{p},'MarkerSize',markersizes{p})
    end
    set(gca,'FontSize',15) 
    xlabel('$t$','FontSize',25,'Interpreter','LaTeX')
    ylabel('$A_t(1,2)$','FontSize',25,'Interpreter','LaTeX') 
    h = legend('True','Uncontrolled','Iteration 1','Iteration 2','Location','Best');
    set(h,'FontSize',16)
    hold off

figure
    hold on
    plot(0:T,Exact.b(:,1),'rx','MarkerSize',20)
    for p = 1:3
        plot(0:T,cSMC{p,1}.b(:,1),'.','Color',colors{p},'MarkerSize',markersizes{p})
    end
    set(gca,'FontSize',15) 
    xlabel('$t$','FontSize',25,'Interpreter','LaTeX')
    ylabel('$b_t(1)$','FontSize',25,'Interpreter','LaTeX') 
    h = legend('True','Uncontrolled','Iteration 1','Iteration 2','Location','Best');
    set(h,'FontSize',16)
    hold off

figure
    hold on
    plot(0:T,Exact.e,'rx','MarkerSize',20)
    for p = 1:3
        plot(0:T,cSMC{p,1}.e,'.','Color',colors{p},'MarkerSize',markersizes{p})
    end
    set(gca,'FontSize',15) 
    xlabel('$t$','FontSize',25,'Interpreter','LaTeX')
    ylabel('$c_t$','FontSize',25,'Interpreter','LaTeX') 
    h = legend('True','Uncontrolled','Iteration 1','Iteration 2','Location','Best');
    set(h,'FontSize',16) 
    hold off
    