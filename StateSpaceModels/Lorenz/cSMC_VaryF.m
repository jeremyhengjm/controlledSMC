% Repeat cSMC as parameter F vary in Lorenz model (d = 8)

d = 2^3;
T = 100;

% Model parameters
sigmasq = 1e-2; % transition noise
sigma = sqrt(sigmasq);
zetasq = 1e-4; % observation noise
zeta = sqrt(zetasq);
interval = 0.1;
sqrtinterval = sqrt(interval);
P = diag(ones(1,d-1),-1);
P(1,d) = 1;
PT = P';
PP = P*P;

% Initial distribution pi_0
Initial = struct();
Initial.Mean = zeros(1,d);
Initial.Cov = sigmasq * eye(d);
Initial.Precision = inv(Initial.Cov);
Initial.MeanPrecision = Initial.Mean * Initial.Precision;
Initial.LogDetCov = log(det(Initial.Cov));
Initial.QuadForm = sum(Initial.MeanPrecision .* Initial.Mean);
Initial.Sample = @(n) mvnrnd(Initial.Mean, Initial.Cov, n); 

% Transition kernel f
Transition = struct();
Transition.G = sigmasq * interval * eye(d); 
Transition.Ginv = inv(Transition.G);
Transition.LogDetG = log(det(Transition.G));

% Load data generated with F = 4.8801
Observation = struct();
load('SimulatedData.mat')
Observation.y = y;
d_y = d-2;
Observation.Std = zeta;
Observation.Precision = 1 / (zeta^2);
    
% Observation density g 
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

% cSMC parameter settings
N = 2^9;
I = 4;
Parameters = struct();
Parameters.Particles = N;
Parameters.Dimension = d;
Parameters.Time = T;    
Parameters.Iterations = I;
Parameters.Approximator = 'quadratic';
Parameters.Adaptive = 'yes';

% Vary F 
ngrid = 40;
gridF = linspace(2.5,8.5,ngrid);
nrep = 100;
cSMC_LogNormConst = zeros(nrep,ngrid);
cSMC_ESS = zeros(nrep,ngrid);
cSMC_Iterations = zeros(nrep,ngrid);

tic
for igrid = 1:ngrid
    igrid
    F = gridF(igrid)
    Lorenz = @(x) ((PT*x' - PP*x').*(P*x') - x' + F)'; 
    Transition.Mean = @(t,x) RungeKutta4(x,interval,Lorenz);

    for irep = 1:nrep
        irep
        cSMC = cSMC_Resample(Parameters, Initial, Transition, Observation, InitialPolicy);
        cSMC_LogNormConst(irep,igrid) = cSMC{end,2}.LogNormConst(end);
        cSMC_ESS(irep,igrid) = mean(cSMC{end,2}.ESS);
        cSMC_Iterations(irep,igrid) = size(cSMC,1) - 1;
    end
    
end
toc

save('cSMC_VaryF.mat','cSMC_LogNormConst','cSMC_ESS','cSMC_Iterations')


