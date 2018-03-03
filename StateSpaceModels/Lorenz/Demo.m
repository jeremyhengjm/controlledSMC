clear 
close all
clc

% load model and data
Model_Lorenz96

%% Initialise with bootstrap particle filter
InitialPolicy = struct();
RepG = zeros(1,d,d);
RepG(1,:,:) = Transition.G;
InitialPolicy.G = repmat(RepG,[T 1 1]);
InitialPolicy.LogDetG = Transition.LogDetG * ones(T,1);
InitialPolicy.A = zeros(T+1,d,d);
InitialPolicy.b = zeros(T+1,d);
InitialPolicy.c = zeros(T+1,1);

%% Initialise with auxiliary particle filter
InitialPolicy = struct();
ProjH = zeros(d_y,d);
ProjH(1:d_y,1:d_y) = eye(d_y);
APF_A = 0.5 * Observation.Precision * (ProjH' * ProjH);
APF_G = inv( Transition.Ginv + 2 * APF_A );
APF_LogDetG = log(det(APF_G));

InitialPolicy.A = zeros(T+1,d,d);
InitialPolicy.b = zeros(T+1,d);
InitialPolicy.c = zeros(T+1,1);
InitialPolicy.G = zeros(T,d,d);
InitialPolicy.LogDetG = zeros(T,1);

InitialPolicy.A(1,:,:) = APF_A;
InitialPolicy.b(1,:) = - Observation.Precision * Observation.y(1,:) * ProjH;
InitialPolicy.c(1) = 0.5 * Observation.Precision * sum(Observation.y(1,:).^2);

for t = 2:(T+1)
    InitialPolicy.A(t,:,:) = APF_A;
    InitialPolicy.b(t,:) = - Observation.Precision * Observation.y(t,:) * ProjH;
    InitialPolicy.c(t) = 0.5 * Observation.Precision * sum(Observation.y(t,:).^2);
    InitialPolicy.G(t-1,:,:) = APF_G;
    InitialPolicy.LogDetG(t-1) = APF_LogDetG;
end

%% Run controlled SMC
N = 2^9;
I = 2;
Parameters = struct();
Parameters.Particles = N;
Parameters.Dimension = d;
Parameters.Time = T;    
Parameters.Iterations = I; % if adapting, this acts as max number of iterations
Parameters.Approximator = 'quadratic';
Parameters.Adaptive = 'no';

tic, cSMC = cSMC_Resample(Parameters, Initial, Transition, Observation, InitialPolicy); toc
   
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
    for i = 1:(I+1)
        IterationNormConsts = cSMC{i,2}.LogNormConst;
        plot(0:T,cSMC{i,2}.LogNormConst,'-*','Color',colors{i})
    end
    set(gca,'FontSize',15) 
    xlabel('$t$','FontSize',25,'Interpreter','LaTeX')
    ylabel('$\log\, Z_t$','FontSize',25,'Interpreter','LaTeX') 
    axis('tight')
    h = InsertLegend(I);
    set(h,'FontSize',16)
            
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
            set(gca,'FontSize',15)
            xlabel('$t$','FontSize',25,'Interpreter','LaTeX')
            ylabel('$\log\,W_t^{\psi^*}(X_{t})$','FontSize',25,'Interpreter','LaTeX')

%% Estimated coefficient A 
figure
    hold on
    plot(0:T,cSMC{2,1}.A(:,1,1),'.','Color',colors{2},'MarkerSize',15)
        for i = 2:(I)
            plot(0:T,cSMC{i+1,1}.A(:,1,1)-cSMC{i,1}.A(:,1,1),'.','Color',colors{i+1},'MarkerSize',15)
        end
        set(gca,'FontSize',15) 
        xlabel('$t$','FontSize',25,'Interpreter','LaTeX')
        ylabel('$A_t^{i}(1,1)$','FontSize',25,'Interpreter','LaTeX') 
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
    ylabel('$b_t^{i}(1)$','FontSize',25,'Interpreter','LaTeX') 
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
        