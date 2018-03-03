function SMC = AIS_Resample(Parameters,Like)
% Annealed importance sampler with Metropolis-adjusted Langevin algorithm

    Grid = Parameters.Grid;
    d = Parameters.Dim;
    N = Parameters.Particles;
    M = Parameters.Steps;
    T = Parameters.TerminalTime;
    P = Parameters.MCMCmoves;
    
    Time = linspace(0,T,M+1);
    StepSize = diff(Time); % 1 x M
    Lambda = linspace(0,1,M+1); 
    DeltaLambda = diff(Lambda); % 1 x M
    
    % Model parameters
    sigmasq = 1.91;
    mu = log(126) - 0.5 * sigmasq;
    beta = 1/33;
    m = 1/d;
    
    % Prior
    Prior = struct();
    Prior.Mean = mu * ones(1,d);
    Prior.Cov = zeros(d,d);
    for i = 1:d
        for j = 1:d
            ind_i = [ floor((i-1)/Grid) + 1; mod(i-1,Grid) + 1 ];
            ind_j = [ floor((j-1)/Grid) + 1; mod(j-1,Grid) + 1 ];
            Prior.Cov(i,j) = sigmasq * exp(- norm(ind_i - ind_j,2) / (Grid * beta) ); 
        end
        
    end
    Prior.Precision = inv(Prior.Cov); 
    Prior.Chol = cholcov(Prior.Cov);
    Prior.MeanPrecision = Prior.Mean * Prior.Precision;
    Prior.Sample = @(n) repmat(Prior.Mean,n,1) + randn(n,d) * Prior.Chol;
    Prior.QuadForm = sum(Prior.MeanPrecision .* Prior.Mean);
    Prior.LogDet = 2 * sum(log(diag(Prior.Chol)));
    Prior.Const = - 0.5 * d * log(2*pi) - 0.5 * Prior.LogDet;
    
    function [f,varargout] = PriorDensity(x,option)
    % Multivariate Gaussian pdf
    % Input: x (N x d)
    %        mu (1 x d)
    %        invSigma (d x d)
    %        const (1 x 1)
    % Output: f (N x 1)
    
        xc = bsxfun(@minus,Prior.Mean,x);
        switch option
            case 'natural'
                f = exp(-0.5 * sum((xc * Prior.Precision) .* xc, 2) + Prior.Const);                
            case 'log'
                f = -0.5 * sum((xc * Prior.Precision) .* xc, 2) + Prior.Const;                
            case 'gradlog'
                f = xc * Prior.Precision;                 
            case 'bothlog'
                xcPrecision = xc * Prior.Precision;
                f = -0.5 * sum(xcPrecision .* xc, 2) + Prior.Const;
                varargout{1} = xcPrecision;     
        end
    end
    Prior.LogDensityAndGrad = @(x) PriorDensity(x,'bothlog');
    
    % Likelihood
    function [l,varargout] = Likelihood(x,option)
        % Likelihood function
        % Input arguments: x (N x d)
        %                  varargin (optional)
        % Output arguments: l (N x 1)
        %                  varargout (gradient, N x d)
        
        mExpX = m * exp(x); % N x d
        mSumExpX = sum(mExpX, 2); % N x 1
        
        switch option
            case 'log'
                l = sum( bsxfun(@times, Like.Obs, x), 2) - mSumExpX; % N x 1
            case 'gradlog'
                l = bsxfun(@minus, Like.Obs, mExpX); % N x d
            case 'bothlog'
                l = sum( bsxfun(@times, Like.Obs, x), 2) - mSumExpX; % N x 1
                varargout{1} = bsxfun(@minus, Like.Obs, mExpX); % N x d

        end
    end
    Like.LogFunc = @(x) Likelihood(x,'log');
    Like.LogFuncAndGrad = @(x) Likelihood(x,'bothlog');
    
    % Target
    Target = struct();
    Target.LogDensity = @(k,x) PriorDensity(x,'log') + Lambda(k) * Likelihood(x,'log');
    Target.GradLogDensity = @(k,x) PriorDensity(x,'gradlog') + Lambda(k) * Likelihood(x,'gradlog');
 
    % Uncontrolled 
    Euler = @(k,x,grad) x + 0.5 * StepSize(k) * grad;
    EulerConst = - 0.5 * d * ( log(2*pi) + log(StepSize) ); % 1 x M
    function f = LogEulerMaruyama(k,EulerMove,x)
    % Log Euler-Maruyama pdf
    % Input:    k (1 x 1)
    %           EulerMove (N x d)
    %           x (N x d)
    % Output: f (N x 1)
    
        xc = x - EulerMove;
        f = - (0.5 / StepSize(k)) * sum(xc .* xc, 2) + EulerConst(k);
    end
    
    function SMC = RunAIS
        % Annealed importance sampler with resampling
            % Output argument: SMC (struct)
            %                 (fields)
            %                   SMC.Trajectory (N x d x M+1)
            %                   SMC.ForwardEuler (N x d x M)
            %                   SMC.LogIncrementalAIS (N x M+1)
            %                   SMC.LogWeights (N x M+1)
            %                   SMC.ESS (1 x M+1)
            %                   SMC.LogNormConst (1 x M+1)

        % Pre-allocate
        SMC = struct();
        SMC.Trajectory = zeros(N,d,M+1);
        SMC.LogIncrementalWeight = zeros(N,M+1);
        SMC.ESS = zeros(1,M+1);
        SMC.LogNormConst = zeros(1,M+1);
        SMC.AvgAcceptProb = zeros(2,M);

        % Initialise
        X = Prior.Sample(N);
        SMC.Trajectory(:,:,1) = X;
        LogRatioNormConst = 0;
        SMC.ESS(1) = N;
        
        for k = 1:M            
            % Weighting            
            LogWeights = DeltaLambda(k) * Like.LogFunc(X);
            SMC.LogIncrementalWeight(:,k+1) = LogWeights;
            MaxLogWeight = max(LogWeights);
            Weights = exp(LogWeights - MaxLogWeight);
            NormalisedWeights = Weights / sum(Weights);
            SMC.ESS(k+1) = 1 / sum(NormalisedWeights.^2);
            LogRatioNormConst = LogRatioNormConst + log(mean(Weights)) + MaxLogWeight;
            SMC.LogNormConst(k+1) = LogRatioNormConst;
            
            % Resample
            Ancestor = Systematic_Resampling(NormalisedWeights);
            X = X(Ancestor,:);
          
            % MCMC move
            AvgAcceptProb = zeros(1,P);
            for p = 1:P
                ForwardEulerMove = Euler(k,X,Target.GradLogDensity(k+1,X)); % N x d 
                Y = mvnrnd(ForwardEulerMove, StepSize(k) * eye(d));
                BackwardEulerMove = Euler(k,Y,Target.GradLogDensity(k+1,Y)); % N x d 
                
                LogMetropolisRatio = Target.LogDensity(k+1,Y) - Target.LogDensity(k+1,X) ... 
                     + LogEulerMaruyama(k,BackwardEulerMove,X) - LogEulerMaruyama(k,ForwardEulerMove,Y); % N x 1
                
                LogAcceptanceProb = LogMetropolisRatio .* (LogMetropolisRatio < 0);
                AcceptLogicals = (log(rand([N,1])) < LogAcceptanceProb);
                X = bsxfun(@times,Y,AcceptLogicals) + bsxfun(@times,X,(1-AcceptLogicals)); 
                AvgAcceptProb(p) = mean(exp(LogAcceptanceProb));
                
            end
            
            SMC.AvgAcceptProb(:,k) = [min(AvgAcceptProb) ; max(AvgAcceptProb)];
            
            % Store
            SMC.Trajectory(:,:,k+1) = X;            
            
        end
        
    end

    SMC = RunAIS;
              
end



