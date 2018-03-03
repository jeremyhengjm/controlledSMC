% Repeat APF and cSMC as observation noise vary in Lorenz model
clear
close all
clc

d = 2^6;
T = 100;

% Model parameters
F = 4.8801; % forcing parameter
sigmasq = 1e-2; % transition noise
sigma = sqrt(sigmasq);
interval = 0.1;
sqrtinterval = sqrt(interval);
P = diag(ones(1,d-1),-1);
P(1,d) = 1;
PT = P';
PP = P*P;
Lorenz = @(x) ((PT*x' - PP*x').*(P*x') - x' + F)';     

% Initial distribution 
Initial = struct();
Initial.Mean = zeros(1,d);
Initial.Cov = sigmasq * eye(d);
Initial.Precision = inv(Initial.Cov);
Initial.MeanPrecision = Initial.Mean * Initial.Precision;
Initial.LogDetCov = log(det(Initial.Cov));
Initial.QuadForm = sum(Initial.MeanPrecision .* Initial.Mean);
Initial.Sample = @(n) mvnrnd(Initial.Mean, Initial.Cov, n); 

% Transition kernel
Transition = struct();
Transition.Mean = @(t,x) RungeKutta4(x,interval,Lorenz);
Transition.G = sigmasq * interval * eye(d); 
Transition.Ginv = inv(Transition.G);
Transition.LogDetG = log(det(Transition.G));
Transition.Sample = @(x,n) mvnrnd(Transition.Mean(1,x), Transition.G, n);

% APF parameter settings
I = 0;
APF_Parameters = struct();
APF_Parameters.Particles = 11468; 
APF_Parameters.Dimension = d;
APF_Parameters.Time = T;    
APF_Parameters.Iterations = I;
APF_Parameters.Approximator = 'quadratic';
APF_Parameters.Adaptive = 'no';

% cSMC parameter settings
I = 1;
cSMC_Parameters = struct();
cSMC_Parameters.Particles = 2^12; 
cSMC_Parameters.Dimension = d;
cSMC_Parameters.Time = T;    
cSMC_Parameters.Iterations = I; 
cSMC_Parameters.Approximator = 'quadratic';
cSMC_Parameters.Adaptive = 'no';

% Loop over zetasq
d_y = d-2;
nzetasq = 3;
gridzetasq = 10.^(-4:-2);
nrep = 100;
APF_LogNormConst = zeros(nrep,nzetasq);
cSMC_LogNormConst = zeros(nrep,nzetasq);
APF_SquaredError = zeros(nrep,nzetasq); 
cSMC_SquaredError = zeros(nrep,nzetasq);
APF_ESS = zeros(nrep,nzetasq);
cSMC_ESS = zeros(nrep,nzetasq);

tic

for izetasq = 2:nzetasq
    % Generate observation
    zetasq = gridzetasq(izetasq);
    zeta = sqrt(zetasq);
    Observation = struct();
    Latentx = zeros(T+1,d);
    Observation.y = zeros(T+1,d_y);
    Latentx(1,:) = Initial.Sample(1);
    Observation.Std = zeta;
    Observation.Precision = 1 / (zetasq);
    Observation.y(1,:) = Latentx(1,1:d_y) + Observation.Std * randn(1,d_y);
    for t = 1:T
        Latentx(t+1,:) = Transition.Sample(Latentx(t,:), 1);
        Observation.y(t+1,:) = Latentx(t+1,1:d_y) + Observation.Std * randn(1,d_y);
    end
    Observation.LogDensity = @(t,x) -0.5 * Observation.Precision * sum( bsxfun(@minus,Observation.y(t,:),x(:,1:d_y)).^2 , 2);

    % Initialise with auxiliary particle filter
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
    
    for irep = 1:nrep
        [izetasq, irep]
        % Run APF
        APF = cSMC_Resample(APF_Parameters, Initial, Transition, Observation, InitialPolicy);
        OutputAPF = APF{1,2};
        APF_LogNormConst(irep,izetasq) = OutputAPF.LogNormConst(end);  
        APF_ESS(irep,izetasq) = mean(OutputAPF.ESS);
        APF_Mean = zeros(T+1,d);
        for s = 1:(T+1)
            CurrentAncestor = OutputAPF.Ancestry(:,s);
            APF_Mean(s,:) = mean(squeeze(OutputAPF.Trajectory(CurrentAncestor,:,s)),1);
        end
        APF_SquaredError(irep,izetasq) = sum( sum( (APF_Mean - Latentx).^2, 2) );
        
        % Run cSMC
        cSMC = cSMC_Resample(cSMC_Parameters, Initial, Transition, Observation, InitialPolicy);
        OutputSMC = cSMC{end,2};
        cSMC_LogNormConst(irep,izetasq) = OutputSMC.LogNormConst(end);
        cSMC_ESS(irep,izetasq) = mean(OutputSMC.ESS);
        TraceAncestor = OutputSMC.Ancestry(:,end);
        cSMC_Mean = zeros(T+1,d);
        cSMC_Mean(end,:) = mean(squeeze(OutputSMC.Trajectory(TraceAncestor,:,end)),1); % compute smoothing mean (1 x d)
        for s = 1:T
            TraceAncestor = OutputSMC.Ancestry(TraceAncestor,s);
            cSMC_Mean(s,:) = mean(squeeze(OutputSMC.Trajectory(TraceAncestor,:,s)),1); % compute smoothing mean (1 x d)
        end
        cSMC_SquaredError(irep,izetasq) = sum( sum( (cSMC_Mean - Latentx).^2, 2) );
        
    end
   
end
toc

save('VaryZetasq_Dim64.mat','APF_LogNormConst','cSMC_LogNormConst','APF_SquaredError','cSMC_SquaredError','APF_ESS','cSMC_ESS')


