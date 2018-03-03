function cSMC = cSMC_Resample(Parameters,Like)
% Iterative SMC sampler
%     clear all
%     close all
%     clc

    Grid = Parameters.Grid;
    d = Parameters.Dim;
    N = Parameters.Particles;
    M = Parameters.Steps;
    T = Parameters.TerminalTime;
    P = Parameters.Iterations;
    Approximator = Parameters.Approximator;
    
    Time = linspace(0,T,M+1);
    StepSize = diff(Time); % 1 x M
    Lambda = linspace(0,1,M+1); 
    DeltaLambda = diff(Lambda); % 1 x M+1 (last element is just a placeholder)
    
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
    Like.LogFuncAndGrad = @(x) Likelihood(x,'bothlog');
 
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
    
    % Controlled
    V0 = @(x,beta) sum((x * beta{1}) .* x,2) ... 
                    + sum(bsxfun(@times,x,beta{2}),2) ... 
                    + beta{3};
    
    function SMC = RunSMC(Twisting)
        % SMC sampler with a given policy
            % Input argument: Optimal (struct)
            %                 (fields)
            %                   Twisting.A (M+1 x d x d)
            %                   Twisting.b (M+1 x d)
            %                   Twisting.d (logical)
            %                   Twisting.e (M+1 x 1)
            %                   Twisting.Theta (M x d x d)
            %                   Twisting.LogDetTheta (M x 1)
            % Output argument: SMC (struct)
            %                 (fields)
            %                   SMC.Trajectory (N x d x M+1)
            %                   SMC.ForwardEuler (N x d x M)
            %                   SMC.ShiftForwardEuler (N x d x M)
            %                   SMC.QuadForm (N x M)
            %                   SMC.LogIncrementalAIS (N x M+1)
            %                   SMC.LogIncrementalWeight (N x M+1)
            %                   SMC.LogWeights (N x M+1)
            %                   SMC.ESS (1 x M+1)
            %                   SMC.LogNormConst (1 x M+1)
            %                   SMC.d (N x M+1)
        InitialPrecision = Prior.Precision + 2 * squeeze(Twisting.A(1,:,:));
        InitialCov = inv(InitialPrecision);
        InitialMean = (Prior.MeanPrecision - Twisting.b(1,:)) * InitialCov;

        % Pre-allocate
        SMC = struct();
        SMC.Trajectory = zeros(N,d,M+1);
        SMC.ShiftForwardEuler = zeros(N,d,M);
        SMC.QuadForm = zeros(N,M);
        SMC.LogIncrementalWeight = zeros(N,M+1);
        SMC.Ancestry = ones(N,M+2);
        SMC.ESS = zeros(1,M+1);
        SMC.LogNormConst = zeros(1,M+1);
        SMC.d = zeros(N,M+1);
        dFunctionNext = zeros(N,1);

        % Initialise
        X = mvnrnd(InitialMean,InitialCov,N);
        SMC.Trajectory(:,:,1) = X;
        [LogPriorCurrent,GradLogPriorCurrent] = Prior.LogDensityAndGrad(X);
        [LogLikeCurrent,GradLogLikeCurrent] = Like.LogFuncAndGrad(X);
        LogDensityCurrent = LogPriorCurrent;
        TransitionGradient = GradLogPriorCurrent + Lambda(2) * GradLogLikeCurrent;
 
        % Compute first transition kernel
        TransitionTheta = squeeze(Twisting.Theta(1,:,:)); % d x d
        TransitionCov = StepSize(1) * TransitionTheta;
        ForwardEulerMove = Euler(1,X,TransitionGradient); % N x d  
        ShiftForwardEuler = bsxfun(@minus,ForwardEulerMove,StepSize(1) * Twisting.b(2,:)); % N x d
        SMC.ShiftForwardEuler(:,:,1) = ShiftForwardEuler;
        TransitionMean = ShiftForwardEuler * TransitionTheta; % N x d
        
        % Compute first weight function
        SMC.V = Twisting.e(1) + 0.5 * Prior.LogDet - 0.5 * log(det(InitialCov)) ... 
                        + 0.5 * (Prior.QuadForm - sum((InitialMean * InitialPrecision) .* InitialMean));
        
        QuadForm = sum((ShiftForwardEuler * TransitionTheta) .* ShiftForwardEuler,2); % N x 1
        SMC.QuadForm(:,1) = QuadForm; 
        
        SMC.d(:,2) = - DeltaLambda(1) * LogLikeCurrent;
%         SMC.d(:,2) = LogDensityCurrent - 0.5 * sum(bsxfun(@plus,TransitionGradient,kappa) .* X,2); % ... 
%                         - 0.125 * StepSize(1) * sum(TransitionGradient .* TransitionGradient,2); % N x 1 
        if (Twisting.d)
            dFunctionCurrent = SMC.d(:,2);            
        else
            dFunctionCurrent = zeros(N,1);
        end
        LogNormalisingFunc = - dFunctionCurrent - Twisting.e(2) + 0.5 * Twisting.LogDetTheta(1) ... 
                             - 0.5 * (1 / StepSize(1)) * ( sum(ForwardEulerMove .* ForwardEulerMove,2) ... 
                                   - QuadForm ); % N x 1
                               
        InitialBeta = {squeeze(Twisting.A(1,:,:)), Twisting.b(1,:), Twisting.e(1)};
        
        LogWeights = - SMC.V + LogNormalisingFunc + V0(X,InitialBeta); 
        
        % Normalize weights, compute ESS and normalizing constant
        MaxLogWeight = max(LogWeights);
        Weights = exp(LogWeights - MaxLogWeight);
        SMC.LogIncrementalWeight(:,1) = LogWeights + SMC.V;
        NormalisedWeights = Weights / sum(Weights);
        SMC.ESS(1) = 1 / sum(NormalisedWeights.^2);
        LogRatioNormConst = log(mean(Weights)) + MaxLogWeight;
        SMC.LogNormConst(1) = LogRatioNormConst;
        
        % Resample
        Ancestor = Systematic_Resampling(NormalisedWeights);
        SMC.Ancestry(:,2) = Ancestor;
        X = X(Ancestor,:);
        ForwardEulerMove = ForwardEulerMove(Ancestor,:);
        LogDensityCurrent = LogDensityCurrent(Ancestor,:);
        dFunctionCurrent = dFunctionCurrent(Ancestor,:);
        TransitionMean = TransitionMean(Ancestor,:);

        
        % Forward controlled simulation
        for k = 1:M
            % Move
            Y = mvnrnd(TransitionMean,TransitionCov);
            [LogPriorNext,GradLogPriorNext] = Prior.LogDensityAndGrad(Y);
            [LogLikeNext,GradLogLikeNext] = Like.LogFuncAndGrad(Y);
            LogDensityNext = LogPriorNext + Lambda(k+1) * LogLikeNext; % log pi_k at Y 
            LogGradNext = GradLogPriorNext + Lambda(k+1) * GradLogLikeNext; % gradient log pi_k at Y
            SMC.Trajectory(:,:,k+1) = Y;
            
            % Uncontrolled weights            
            BackwardEulerMove = Euler(k,Y,LogGradNext); 
            LogIncrementalAIS = LogEulerMaruyama(k,BackwardEulerMove,X) - LogEulerMaruyama(k,ForwardEulerMove,Y) ... 
                              + LogDensityNext - LogDensityCurrent;                                                   
            
            if (k < M)
                % Compute next transition kernel
                TransitionTheta = squeeze(Twisting.Theta(k+1,:,:)); % d x d
                TransitionCov = StepSize(k+1) * TransitionTheta;
                TransitionGradient = LogGradNext + DeltaLambda(k+1) * GradLogLikeNext; % gradient log pi_{k+1} at Y
                ForwardEulerMove = Euler(k+1,Y,TransitionGradient); % N x d  
                ShiftForwardEuler = bsxfun(@minus,ForwardEulerMove,StepSize(k+1) * Twisting.b(k+2,:)); % N x d
                SMC.ShiftForwardEuler(:,:,k+1) = ShiftForwardEuler;
                TransitionMean = ShiftForwardEuler * TransitionTheta; % N x d

                % Compute weight function
                QuadForm = sum((ShiftForwardEuler * TransitionTheta) .* ShiftForwardEuler,2); % N x 1
                SMC.QuadForm(:,k+1) = QuadForm; 

                SMC.d(:,k+2) = -DeltaLambda(k+1) * LogLikeNext;
%                 SMC.d(:,k+2) = LogDensityNext - 0.5 * sum(bsxfun(@plus,TransitionGradient,kappa) .* Y,2); % ... 
%                                 - 0.125 * StepSize(k+1) * sum(TransitionGradient .* TransitionGradient,2); % N x 1  
                if (Twisting.d)
                    dFunctionNext = SMC.d(:,k+2);
                else
                    dFunctionNext = zeros(N,1);
                end                
                LogNormalisingFunc = - dFunctionNext - Twisting.e(k+2) + 0.5 * Twisting.LogDetTheta(k+1) ... 
                                     - 0.5 * (1 / StepSize(k+1)) * ( sum(ForwardEulerMove .* ForwardEulerMove,2) ... 
                                           - QuadForm ); % N x 1
            else
                LogNormalisingFunc = zeros(N,1);
            end

            TransitionBeta = {squeeze(Twisting.A(k+1,:,:)), Twisting.b(k+1,:), Twisting.e(k+1)};
                            
            currentV = V0(Y,TransitionBeta) + dFunctionCurrent; % N x 1
            
            LogWeights = LogIncrementalAIS + LogNormalisingFunc + currentV; % N x 1
            SMC.LogIncrementalWeight(:,k+1) = LogWeights;
                        
            % Normalize weights, compute ESS and normalizing constant
            MaxLogWeight = max(LogWeights);
            Weights = exp(LogWeights - MaxLogWeight);
            NormalisedWeights = Weights / sum(Weights);
            SMC.ESS(k+1) = 1 / sum(NormalisedWeights.^2);
            LogRatioNormConst = LogRatioNormConst + log(mean(Weights)) + MaxLogWeight;
            SMC.LogNormConst(k+1) = LogRatioNormConst;

            % Resampling
            Ancestor = Systematic_Resampling(NormalisedWeights);
            SMC.Ancestry(:,k+2) = Ancestor;
            X = Y(Ancestor,:);
            ForwardEulerMove = ForwardEulerMove(Ancestor,:);
            LogDensityCurrent = LogDensityNext(Ancestor,:);
            dFunctionCurrent = dFunctionNext(Ancestor,:);
            TransitionMean = TransitionMean(Ancestor,:);
            
        end
        
    end

    % Pre-compute
    switch Approximator
        case 'quadratic'
            refA = d * (d + 1) / 2 + 1;
            IndexA = (d + 2):refA;
            IndexAA = (refA+1):(refA+d); 
            
        case 'purequadratic'
            IndexAA = (d + 2):(2 * d + 1); 
            
    end
    Indexb = 2:(d+1);
    Indexe = 1;
    
    function NewPolicy = ApproxDP(SMC,Policy)
        % Approximate dynamic programming 
            % Input argument: SMC (struct)
            %                 (fields)
            %                   SMC.Trajectory (N x d x M+1)
            %                   SMC.ShiftForwardEuler (N x d x M)
            %                   SMC.QuadForm (N x M)
            %                   SMC.LogIncrementalWeight (N x M+1)
            %                   SMC.d (N x M+1)
            %                   SMC.Ancestry (N x M+2)
            %                 Policy (struct)
            %                 (fields)
            %                   Policy.A (M+1 x d x d)
            %                   Policy.b (M+1 x d)
            %                   Policy.d (logical)
            %                   Policy.e (M+1 x 1)
            %                   Policy.Theta (M x d x d)
            %                   Policy.LogDetTheta (M x 1)
            %                 First (string)
            %                 WeightedLS (string)
            % Output argument: SMC (struct)
            %                 NewPolicy (struct)
            %                 (fields)
            %                   NewPolicy.A (M+1 x d x d)
            %                   NewPolicy.b (M+1 x d)
            %                   NewPolicy.d (logical)
            %                   NewPolicy.e (M+1 x 1)
            %                   NewPolicy.Theta (M x d x d)
            %                   NewPolicy.LogDetTheta (M x 1)
            
  
        % Pre-allocate
        NewPolicy = struct();
        NewPolicy.A = zeros(M+1,d,d);
        NewPolicy.b = zeros(M+1,d);
        NewPolicy.e = zeros(M+1,1);
        NewPolicy.Theta = zeros(M,d,d);
        NewPolicy.LogDetTheta = zeros(M,1);
        NewPolicy.d = true;
        
        % Initialise
        LogNormalisingFunc = zeros(N,1);
        if (Policy.d)
            dFunction = zeros(N,M+1);
        else
            dFunction = SMC.d;
        end

        for k = (M+1):-1:1
            % Compute value at particle locations
            Value = - SMC.LogIncrementalWeight(:,k) - LogNormalisingFunc - dFunction(SMC.Ancestry(:,k),k); % N x 1
            
            DesignMatrix = x2fx(squeeze(SMC.Trajectory(:,:,k)),Approximator);
            if (k == 1)
                BetaHat = (DesignMatrix' * DesignMatrix) \ (DesignMatrix' * Value); % (d^2/2 + 3d/2 + 1) x 1
                
                if ( d == 1 || strcmp(Approximator,'purequadratic') )
                    NewPolicy.A(k,:,:) = squeeze(Policy.A(k,:,:)) + diag(BetaHat(IndexAA));
                else
                    NewPolicy.A(k,:,:) = squeeze(Policy.A(k,:,:)) + diag(BetaHat(IndexAA)) + 0.5 * squareform(BetaHat(IndexA));
                end
                NewPolicy.b(k,:) = Policy.b(k,:) + BetaHat(Indexb)';   
                NewPolicy.e(k) = Policy.e(k) + BetaHat(Indexe);
                
            else
                % Fit value function (least squares estimation)
                BetaHat = (DesignMatrix' * DesignMatrix) \ (DesignMatrix' * Value); % (d^2/2 + 3d/2 + 1) x 1
                
                % Store parameters
                if ( d == 1 || strcmp(Approximator,'purequadratic') )
                    AHat = diag(BetaHat(IndexAA));
                else
                    AHat = diag(BetaHat(IndexAA)) + 0.5 * squareform(BetaHat(IndexA));
                end
                bHat = BetaHat(Indexb)';
                eHat = BetaHat(Indexe);
                
                UpdatedA = squeeze(Policy.A(k,:,:)) + AHat;
                Theta = inv( eye(d) + 2 * StepSize(k-1) * UpdatedA ); % updated
                LogDetTheta = log(det(Theta));

                NewPolicy.A(k,:,:) = UpdatedA;
                NewPolicy.b(k,:) = Policy.b(k,:) + bHat;
                NewPolicy.e(k) = Policy.e(k) + eHat; 
                NewPolicy.Theta(k-1,:,:) = Theta;
                NewPolicy.LogDetTheta(k-1) = LogDetTheta;

                % Evaluate next normalising function   
                NewShiftForwardEuler = bsxfun(@minus,squeeze(SMC.ShiftForwardEuler(:,:,k-1)),StepSize(k-1) * bHat); % N x d
                LogNormalisingFunc = - dFunction(:,k) - eHat + 0.5 * LogDetTheta - 0.5 * Policy.LogDetTheta(k-1) ... 
                                     - 0.5 * (1 / StepSize(k-1)) * ( SMC.QuadForm(:,k-1) ... 
                                           - sum((NewShiftForwardEuler * Theta) .* NewShiftForwardEuler,2) ); % N x 1
                
            end
            
        end        
             
    end
         
    % Pre-allocate
    cSMC = cell(P+1,2); % store value functions (left) and SMC output (right)
    
    % Initialise with AIS
    CurrentPolicy = struct();
    Id = zeros(1,d,d);
    Id(1,:,:) = eye(d);
    CurrentPolicy.Theta = repmat(Id,[M 1 1]);
    CurrentPolicy.LogDetTheta = zeros(M,1);
    CurrentPolicy.A = zeros(M+1,d,d);
    CurrentPolicy.b = zeros(M+1,d);
    CurrentPolicy.d = false;
    CurrentPolicy.e = zeros(M+1,1);
    cSMC{1,1} = CurrentPolicy;
    OutputSMC = RunSMC(CurrentPolicy); 
    cSMC{1,2} = OutputSMC;

    % Iterate 
    for p = 1:P
        UpdatedPolicy = ApproxDP(OutputSMC,CurrentPolicy);         
        cSMC{p+1,1} = UpdatedPolicy;
        CurrentPolicy = UpdatedPolicy;        
        OutputSMC = RunSMC(CurrentPolicy); 
        cSMC{p+1,2} = OutputSMC;
    end   
        
end



