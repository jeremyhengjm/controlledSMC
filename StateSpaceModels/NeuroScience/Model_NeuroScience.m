% Demba's neuroscience example

d = 1;

% Model parameters
alpha = 0.99;
sigmasq = 0.11;

% Initial distribution 
Initial = struct();
Initial.Mean = 0;
Initial.Cov = 1;
Initial.Precision = 1;
Initial.MeanPrecision = Initial.Mean * Initial.Precision;
Initial.LogDetCov = log(Initial.Cov);
Initial.QuadForm = Initial.MeanPrecision * Initial.Mean;
Initial.Sample = @(n) mvnrnd(Initial.Mean, Initial.Cov, n); 

% Transition kernel 
Transition = struct();
Transition.F = alpha;
Transition.Mean = @(t,x) x * Transition.F;
Transition.G = sigmasq;
Transition.Ginv = 1 / sigmasq;
Transition.LogDetG = log(Transition.G);
Transition.Sample = @(x,n) mvnrnd(Transition.Mean(1,x), Transition.G, n);
Transition.Const = -0.5 * log(2 * pi * sigmasq);
Transition.LogDensity = @(x,y) - 0.5 * (y - alpha * x).^2 / sigmasq; % + Transition.Const;

% Observation
expit = @(x) exp(x) ./ (1 + exp(x));
ntrials = 50;
Observation = struct();    

% Load data
load('Data_NeuroScience.mat')
% load('SimulatedData_NeuroScience.mat')
Observation.y = y;
T = length(y) - 1;
Observation.Const = zeros(T+1,1);
for t = 1:(T+1)
    Observation.Const(t) = log( nchoosek(ntrials, Observation.y(t)) );
end
figure
    hold on
    plot(1:(T+1), Observation.y,'b-')
    axis('tight')
    set(gca,'FontSize',15) % tick size
    xlabel('$t$','Interpreter','LaTeX','FontSize',25)
    ylabel('$y_t$','Interpreter','LaTeX','FontSize',25)
    
% Observation density g 
Observation.LogDensity = @(t,x) Observation.Const(t) + Observation.y(t) * log( expit(x) ) ... 
                                + ( ntrials - Observation.y(t) ) * log( 1 - expit(x) );
