% Lorenz model

d = 2^3; % dimension 
T = 100; % no. of time steps

% Model parameters
F = 4.8801; % forcing parameter
sigmasq = 1e-2; % transition noise
sigma = sqrt(sigmasq);
zetasq = 1e-2; % observation noise
zeta = sqrt(zetasq);
interval = 0.1;
sqrtinterval = sqrt(interval);

P = diag(ones(1,d-1),-1);
P(1,d) = 1;
PT = P';
PP = P*P;

Lorenz = @(x) ((PT*x' - PP*x').*(P*x') - x' + F)'; % drift

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

% Load generated observation for reproducibility
Observation = struct();
load('SimulatedData_SmallNoise.mat')
Observation.y = y;
d_y = d-2;
Observation.Std = zeta;
Observation.Precision = 1 / (zeta^2);

% Generate obs
% Observation = struct();
% Latentx = zeros(T+1,d);
% d_y = d-2;
% Observation.y = zeros(T+1,d_y);
% Latentx(1,:) = Initial.Sample(1);
% Observation.Std = zeta;
% Observation.Precision = 1 / (zetasq);
% Observation.y(1,:) = Latentx(1,1:d_y) + Observation.Std * randn(1,d_y);
% for t = 1:T
%     Latentx(t+1,:) = Transition.Sample(Latentx(t,:), 1);
%     Observation.y(t+1,:) = Latentx(t+1,1:d_y) + Observation.Std * randn(1,d_y);
% end
% y = Observation.y;
% save('SimulatedData_new.mat','Latentx','y')

figure
    for i = 1:4
        subplot(2,2,i)
        hold on
        plot(1:(T+1), Observation.y(:,i),'b-*')
        plot(1:(T+1), Latentx(:,i),'r-*')        
        axis('tight')
        set(gca,'FontSize',15) % tick size
        xlabel('$t$','Interpreter','LaTeX','FontSize',25)
        h = legend('$x_t$','$y_t$');
        set(h,'Interpreter','LaTeX','FontSize',25)
    end
    
% Observation density g 
% Observation.Const = - 0.5 * d_y * ( log(2*pi) + log(zeta^2) );
Observation.LogDensity = @(t,x) -0.5 * Observation.Precision * sum( bsxfun(@minus,Observation.y(t,:),x(:,1:d_y)).^2 , 2);% + Observation.Const;
