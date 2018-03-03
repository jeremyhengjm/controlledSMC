function cSMC = LQG(Parameters)
% Optimal LQG controlled SMC 

    d = Parameters.Dim;
    N = Parameters.Particles;
    M = Parameters.Steps;
    T = Parameters.TerminalTime;
    Time = linspace(0,T,M+1);
    StepSize = diff(Time); % 1 x M
    Lambda = linspace(0,1,M+1); 
    DeltaLambda = diff(Lambda); % 1 x M+1 (last element is just a placeholder)
    
    % Prior
    Prior = struct();
    Prior.Mean = 0 * ones(1,d);
    Prior.Cov = eye(d);
    Prior.Chol = cholcov(Prior.Cov);
    Prior.Precision = inv(Prior.Cov);
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
    
    % Likelihood (Model: y = G(x) + noise)
    Like = struct();
    Like.G = @(x) x;
    Like.Obs = Parameters.Obs;
    Like.Rho = Parameters.Rho;
    Like.Var = 1 * ones(1,d); 
    Like.Cov = Like.Rho * ones(d,d) + diag(Like.Var - Like.Rho); % size: d x d
    Like.Precision = inv(Like.Cov);
    Like.QuadForm = sum((Like.Obs * Like.Precision) .* Like.Obs);
    function [l,varargout] = Likelihood(x,option)
        % Likelihood function
        
        xc = bsxfun(@minus,Like.Obs,Like.G(x)); % N x d
        switch option
            case 'natural'
                l = exp(-0.5 * sum((xc * Like.Precision) .* xc, 2));                
            case 'log'
                l = -0.5 * sum((xc * Like.Precision) .* xc, 2);                
            case 'gradlog'
                l = xc * Like.Precision;                
            case 'bothlog'
                xcPrecision = xc * Like.Precision;
                l = -0.5 * sum(xcPrecision .* xc, 2); % log-scale
                varargout{1} = xcPrecision;
        end
    end
    Like.LogFuncAndGrad = @(x) Likelihood(x,'bothlog');
    
    % True target
    True = struct();
    True.Precision = @(k) Prior.Precision + Lambda(k) * Like.Precision;
    True.Cov = @(k) inv(True.Precision(k));
    True.Mean = @(k) (True.Cov(k) * (Prior.Precision * (Prior.Mean') + Lambda(k) * Like.Precision * (Like.Obs')))'; 
    True.Density = @(k,x) mvnpdf(x,True.Mean(k),True.Cov(k)); 
    True.LogDensity = @(k,x) logmvnpdf(x,True.Mean(k),True.Cov(k));
    True.LogNormConst = @(k) -0.5 * Prior.LogDet + 0.5 * log(det(True.Cov(k))) ... 
        - 0.5 * Prior.QuadForm ... 
        - 0.5 * Lambda(k) * Like.QuadForm ... 
        + 0.5 * sum((True.Mean(k) * True.Precision(k)) .* True.Mean(k));
    True.v = @(k) Prior.Mean * Prior.Precision + Lambda(k) * Like.Obs * Like.Precision; % 1 x d
    True.G = @(k) eye(d) - 0.5 * StepSize(k-1) * True.Precision(k); % d x d
    True.h = @(k) 0.5 * StepSize(k-1) * True.v(k); % 1 x d 
 
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
                
    V = @(x,y,beta) sum((y * beta{1}) .* y,2) + sum(bsxfun(@times,y,beta{2}),2) ... 
                     + sum((x * beta{3}) .* x,2) + sum(bsxfun(@times,x,beta{4}),2) ...           
                     + beta{5}; % N x 1
    
    function SMC = RunSMC(Twisting)
        % SMC sampler with a given policy
            % Input argument: Optimal (struct)
            %                 (fields)
            %                   Twisting.A (M+1 x d x d)
            %                   Twisting.b (M+1 x d)
            %                   Twisting.C (M+1 x d x d)
            %                   Twisting.d (M+1 x d)
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
        
        InitialPrecision = Prior.Precision + 2 * squeeze(Twisting.A(1,:,:));
        InitialCov = inv(InitialPrecision);
        InitialMean = (Prior.MeanPrecision - Twisting.b(1,:)) * InitialCov;

        % Pre-allocate
        SMC = struct();
        SMC.Trajectory = zeros(N,d,M+1);
        SMC.ForwardEuler = zeros(N,d,M);
        SMC.ShiftForwardEuler = zeros(N,d,M);
        SMC.QuadForm = zeros(N,M);
        SMC.LogIncrementalAIS = zeros(N,M+1);
        SMC.LogIncrementalWeight = zeros(N,M+1);
        SMC.ESS = zeros(1,M+1);
        SMC.LogNormConst = zeros(1,M+1);

        % Initialise
        X = mvnrnd(InitialMean,InitialCov,N);
        SMC.Trajectory(:,:,1) = X;
        [LogPriorCurrent,GradLogPriorCurrent] = Prior.LogDensityAndGrad(X);
        GradLogLikeCurrent = Likelihood(X,'gradlog');
        LogDensityCurrent = LogPriorCurrent;
        TransitionGradient = GradLogPriorCurrent + Lambda(2) * GradLogLikeCurrent;
 
        % Compute first transition kernel
        TransitionTheta = squeeze(Twisting.Theta(1,:,:)); % d x d
        TransitionCov = StepSize(1) * TransitionTheta;
        ForwardEulerMove = Euler(1,X,TransitionGradient); % N x d  
        SMC.ForwardEuler(:,:,1) = ForwardEulerMove;
        ShiftForwardEuler = bsxfun(@minus,ForwardEulerMove,StepSize(1) * Twisting.b(2,:)); % N x d
        SMC.ShiftForwardEuler(:,:,1) = ShiftForwardEuler;
        TransitionMean = ShiftForwardEuler * TransitionTheta; % N x d
        
        % Compute first weight function
        SMC.V = Twisting.e(1) + 0.5 * Prior.LogDet - 0.5 * log(det(InitialCov)) ... 
                        + 0.5 * (Prior.QuadForm - sum((InitialMean * InitialPrecision) .* InitialMean));
        
        QuadForm = sum((ShiftForwardEuler * TransitionTheta) .* ShiftForwardEuler,2); % N x 1
        SMC.QuadForm(:,1) = QuadForm; 
        LogNormalisingFunc = - sum((X * squeeze(Twisting.C(2,:,:))) .* X,2) - sum(bsxfun(@times,X,Twisting.d(2,:)),2) ... 
                             - Twisting.e(2) + 0.5 * Twisting.LogDetTheta(1) ... 
                             - 0.5 * (1 / StepSize(1)) * ( sum(ForwardEulerMove .* ForwardEulerMove,2) ... 
                                   - QuadForm ); % N x 1
                               
        InitialBeta = {squeeze(Twisting.A(1,:,:)), Twisting.b(1,:), Twisting.e(1)};
        
        LogWeights = - SMC.V + LogNormalisingFunc + V0(X,InitialBeta); 
        
        % Normalize weights, compute ESS and normalizing constant
        MaxLogWeight = max(LogWeights);
        Weights = exp(LogWeights - MaxLogWeight);
        SMC.LogIncrementalWeight(:,1) = LogWeights;
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
        TransitionMean = TransitionMean(Ancestor,:);
        
        % Forward controlled simulation
        for k = 1:M
            % Move
            Y = mvnrnd(TransitionMean,TransitionCov);
            [LogPriorNext,GradLogPriorNext] = Prior.LogDensityAndGrad(Y);
            [LogLikeNext,GradLogLikeNext] = Like.LogFuncAndGrad(Y);
            LogDensityNext = LogPriorNext + Lambda(k+1) * LogLikeNext; % log pi_k at Y 
            LogGradNext = GradLogPriorNext + Lambda(k+1) * GradLogLikeNext; % gradient log pi_k at Y
            
            % Uncontrolled weights            
            BackwardEulerMove = Euler(k,Y,LogGradNext); 
            LogIncrementalAIS = LogEulerMaruyama(k,BackwardEulerMove,X) - LogEulerMaruyama(k,ForwardEulerMove,Y) ... 
                              + LogDensityNext - LogDensityCurrent;                                                   
            SMC.LogIncrementalAIS(:,k+1) = LogIncrementalAIS;
            
            if (k < M)
                % Compute next transition kernel
                TransitionTheta = squeeze(Twisting.Theta(k+1,:,:)); % d x d
                TransitionCov = StepSize(k+1) * TransitionTheta;
                TransitionGradient = LogGradNext + DeltaLambda(k+1) * GradLogLikeNext; % gradient log pi_{k+1} at Y
                ForwardEulerMove = Euler(k+1,Y,TransitionGradient); % N x d  
                SMC.ForwardEuler(:,:,k+1) = ForwardEulerMove;
                ShiftForwardEuler = bsxfun(@minus,ForwardEulerMove,StepSize(k+1) * Twisting.b(k+2,:)); % N x d
                SMC.ShiftForwardEuler(:,:,k+1) = ShiftForwardEuler;
                TransitionMean = ShiftForwardEuler * TransitionTheta; % N x d

                % Compute weight function
                QuadForm = sum((ShiftForwardEuler * TransitionTheta) .* ShiftForwardEuler,2); % N x 1
                SMC.QuadForm(:,k+1) = QuadForm; 

                LogNormalisingFunc = - sum((Y * squeeze(Twisting.C(k+2,:,:))) .* Y,2) - sum(bsxfun(@times,Y,Twisting.d(k+2,:)),2) ... 
                                     - Twisting.e(k+2) + 0.5 * Twisting.LogDetTheta(k+1) ... 
                                     - 0.5 * (1 / StepSize(k+1)) * ( sum(ForwardEulerMove .* ForwardEulerMove,2) ... 
                                           - QuadForm ); % N x 1
            else
                LogNormalisingFunc = zeros(N,1);
            end

            TransitionBeta = {squeeze(Twisting.A(k+1,:,:)), Twisting.b(k+1,:), squeeze(Twisting.C(k+1,:,:)), ... 
                                Twisting.d(k+1,:), Twisting.e(k+1)};
                            
            currentV = V(X,Y,TransitionBeta); % N x 1
            
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
            TransitionMean = TransitionMean(Ancestor,:);
            
        end
        
    end

    function Exact = ComputeLQG
    % Output argument: Exact (struct)
    %                  (fields)
    %                   Exact.A (M+1 x d x d)
    %                   Exact.b (M+1 x d)
    %                   Exact.C (M+1 x d x d)
    %                   Exact.d (M+1 x d)
    %                   Exact.e (M+1 x 1)
    %                   Exact.Theta (M x d x d)
    %                   Exact.LogDetTheta (M x 1)
    
    % Pre-allocate
        Exact.A = zeros(M+1,d,d);
        Exact.b = zeros(M+1,d);
        Exact.C = zeros(M+1,d,d);
        Exact.d = zeros(M+1,d);
        Exact.e = zeros(M+1,1);
        Exact.Theta = zeros(M,d,d);
        Exact.LogDetTheta = zeros(M,1);

    % Initialise
        Exact.A(end,:,:) = 0.125 * StepSize(M) * True.Precision(M+1) * True.Precision(M+1); 
        Exact.b(end,:) = - 0.25 * StepSize(M) * True.v(M+1) * True.Precision(M+1);
        Exact.C(end,:,:) = - squeeze(Exact.A(end,:,:)) + 0.5 * DeltaLambda(M) * Like.Precision;
        Exact.d(end,:) = - DeltaLambda(M) * Like.Obs * Like.Precision - Exact.b(end,:);
        Exact.e(end) = 0.5 * DeltaLambda(M) * Like.QuadForm;

        % Backward recursion
        for k = M:-1:1
            % Precompute
            if (k == 1)
                currentA = zeros(d,d);
                currentb = zeros(1,d);
                currentC = zeros(d,d);
                currentd = zeros(1,d);
                currente = 0;
            else
                currentPrecision = True.Precision(k);
                currentA = 0.125 * StepSize(k-1) * currentPrecision * currentPrecision; % d x d
                currentb = - 0.25 * StepSize(k-1) * True.v(k) * currentPrecision; % 1 x d
                currentC = - currentA + 0.5 * DeltaLambda(k-1) * Like.Precision;
                currentd = - DeltaLambda(k-1) * Like.Obs * Like.Precision - currentb;
                currente = 0.5 * DeltaLambda(k-1) * Like.QuadForm;
            end 

            % Assign
            futureTheta = inv( eye(d) + 2 * StepSize(k) * squeeze(Exact.A(k+1,:,:)) );
            Exact.Theta(k,:,:) = futureTheta;
            Exact.LogDetTheta(k) = log(det(futureTheta));
            futurev = True.v(k+1);
            futureG = True.G(k+1);
            futureh = True.h(k+1);
            futurehDb = futureh - StepSize(k) * Exact.b(k+1,:); % 1 x d
            Exact.A(k,:,:) = currentA + squeeze(Exact.C(k+1,:,:)) ... 
                                + 0.5 * (1 / StepSize(k)) * futureG * (eye(d) - futureTheta) * futureG;
            Exact.b(k,:) = currentb + Exact.d(k+1,:) ... 
                                + 0.5 * futurev * (eye(d) - futureTheta) * futureG ... 
                                + Exact.b(k+1,:) * futureTheta * futureG;
            Exact.C(k,:,:) = currentC;
            Exact.d(k,:) = currentd;
            Exact.e(k) = currente + Exact.e(k+1) - 0.5 * Exact.LogDetTheta(k) ... 
                            + 0.5 * (1 / StepSize(k)) * ( sum(futureh .* futureh) ... 
                            - sum((futurehDb * futureTheta) .* futurehDb) ); 
        end     
    end 
         
    % Pre-allocate
    cSMC = cell(2,2); % store value functions (left) and SMC output (right)
    
    % Initialise with AIS
    CurrentPolicy = struct();
    Id = zeros(1,d,d);
    Id(1,:,:) = eye(d);
    CurrentPolicy.Theta = repmat(Id,[M 1 1]);
    CurrentPolicy.LogDetTheta = zeros(M,1);
    CurrentPolicy.A = zeros(M+1,d,d);
    CurrentPolicy.b = zeros(M+1,d);
    CurrentPolicy.C = zeros(M+1,d,d);
    CurrentPolicy.d = zeros(M+1,d);
    CurrentPolicy.e = zeros(M+1,1);
    cSMC{1,1} = CurrentPolicy;
    OutputSMC = RunSMC(CurrentPolicy); 
    OutputSMC.Diff = OutputSMC.LogNormConst(M+1) - True.LogNormConst(M+1);
    cSMC{1,2} = OutputSMC;
    
    % Run LQG and optimally controlled SMC
    Exact = ComputeLQG;
    Exact.LogNormConst = zeros(1,M+1);
    for m = 1:(M+1)
        Exact.LogNormConst(m) = True.LogNormConst(m);
    end
    cSMC{2,1} = Exact;
    OutputSMC = RunSMC(Exact);
    cSMC{2,2} = OutputSMC;
            
end



