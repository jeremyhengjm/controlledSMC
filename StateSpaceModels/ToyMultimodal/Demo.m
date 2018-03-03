clear 
close all
clc

Model_ToyMultimodal

%% Run controlled SMC
N = 2^9;
M = 2^4;
I = 1;
Parameters = struct();
Parameters.Particles = N;
Parameters.Knots = M;
Parameters.Time = T;    
Parameters.Iterations = I;
Parameters.KnotsType = 'particles';
Parameters.BandwidthFactor = 0.3;

tic, cSMC = cSMC_Resample(Parameters, Initial, Transition, Observation); toc

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
    set(h,'FontSize',15)
%% Plot standard deviation of log weights
figure
    hold on
    for i = 1:(I+1)
        plot(0:T,std(cSMC{i,2}.LogIncrementalWeight),'-*','Color',colors{i})
    end
    set(gca,'FontSize',15)
    xlabel('$t$','FontSize',25,'Interpreter','LaTeX')
    ylabel('Std of log weights','FontSize',25,'Interpreter','LaTeX') 
    h = InsertLegend(I);
    set(h,'FontSize',15)
    
%% Normalising constant estimation
figure
    hold on
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
        set(gca,'FontSize',15)
        xlabel('$t$','FontSize',25,'Interpreter','LaTeX')
        ylabel('$\log\,W_t^{\psi}(X_{t})$','FontSize',25,'Interpreter','LaTeX')
    
%% Plot policy functions
npoints = 1000;
grid = linspace(-50,50,npoints)';
g = @(t,x) exp(Observation.LogDensity(t,x));
RBF = @(beta,distancesq) exp(- beta * distancesq);

figure
    plot(grid, g(T+1,grid), 'r-','LineWidth',1)
    hold on
    plot(grid, EvaluatePolicy(T+1,grid,cSMC{2,1}), 'b-')
    set(gca,'FontSize',13) 
    xlabel('$x_T$','FontSize',25,'Interpreter','LaTeX')
    ylabel('$\psi_T$','FontSize',25,'Interpreter','LaTeX') 
    h = legend('Optimal','Approximation');
    set(h,'FontSize',16)

grid = linspace(-20,20,npoints)';
f = @(t,x,y) normpdf(y, Transition.Mean(t,x), sigmav);
    
NormalisingFunc = zeros(1,npoints);
for ipoint = 1:npoints
    NormalisingFunc(ipoint) = trapz(grid, f(T, grid(ipoint), grid) .* g(T+1, grid) );
end

figure
    plot(grid, g(T, grid)' .* NormalisingFunc, 'r-','LineWidth',1)
    hold on
    plot(grid, EvaluatePolicy(T,grid,cSMC{2,1}), 'b-')
    set(gca,'FontSize',15) % 
    xlabel('$x_{T-1}$','FontSize',25,'Interpreter','LaTeX')
    ylabel('$\psi_{T-1}$','FontSize',25,'Interpreter','LaTeX') 
    h = legend('Optimal','Approximation');
    set(h,'FontSize',16)
