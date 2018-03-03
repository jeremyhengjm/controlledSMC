% Toy multimodal example from Gordon et al (1993) and Kitagawa (1996)

d = 1;
T = 100;

% Model parameters
sigmavsq = 10;
sigmav = sqrt(sigmavsq);
sigmawsq = 1;
sigmaw = sqrt(sigmawsq);

% Initial distribution 
Initial = struct();
Initial.Mean = 0;
Initial.Var = 5;
Initial.Std = sqrt(Initial.Var);
Initial.Precision = 1 / Initial.Var;
Initial.MeanPrecision = Initial.Mean * Initial.Precision;
Initial.LogDetCov = log(Initial.Var);
Initial.QuadForm = Initial.MeanPrecision * Initial.Mean;
Initial.Sample = @(n) Initial.Mean + Initial.Std * randn(n,1);

% Transition kernel 
Transition = struct();
Transition.Mean = @(t,x) 0.5 * x + 25 * x ./ (1 + x.^2) + 8 * cos(1.2 * t);
Transition.G = sigmavsq;
Transition.Ginv = 1 / sigmavsq;
Transition.LogG = log(Transition.G);
Transition.Sample = @(t,x,n) Transition.Mean(t,x)  + sigmav * randn(n,1);

% Observation
Observation = struct();
Observation.Func = @(x) 0.05 * x.^2;
Observation.Std = sigmaw;
Observation.Precision = 1 / Observation.Std^2;

% Load obs
% load('SimulatedData_SigmawSq_1.mat')
% Observation.y = y;

% Generate obs
Latentx = zeros(T+1,1);
Observation.y = zeros(T+1,1);
Latentx(1) = Initial.Sample(1);
Observation.y(1) = Observation.Func(Latentx(1)) + Observation.Std * randn(1);
for t = 1:T
    Latentx(t+1) = Transition.Sample(t, Latentx(t,:), 1);
    Observation.y(t+1) = Observation.Func(Latentx(t+1)) + Observation.Std * randn(1);
end
 
figure
    hold on
    plot(1:(T+1), Observation.y,'b-*')
    axis('tight')
    xlabel('$t$','Interpreter','LaTeX','FontSize',25)
    ylabel('$y_t$','Interpreter','LaTeX','FontSize',25)
    
% Observation density g 
Observation.Const = - 0.5 * log(2 * pi * sigmawsq);
Observation.LogDensity = @(t,x) -0.5 * Observation.Precision * ( Observation.y(t) - Observation.Func(x) ).^2 + Observation.Const;
