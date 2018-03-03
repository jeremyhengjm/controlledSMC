function cSMC = FFBS(Parameters, Initial, Transition, Observation)
% Forward filtering backward smoothing
    % Input arguments: Parameters (struct)
    %                  Initial (struct)
    %                  Transition (struct)
    %                  Observation (struct)
    % Output arguments: cSMC (struct)
    %                   (fields)
    %                   Trajectory (N x T+1)
    %                   SmoothingWeights (N x T+1) 

    N = Parameters.Particles;
    d = Parameters.Dimension;
    T = Parameters.Time;
    
    % input: x (N x d), A (d x d), output: N x 1
    QuadraticForm = @(x,A) sum((x * A) .* x,2); 
    % input: x (N x d), beta{1} (d x d), beta{2} (1 x d), beta{3} (1 x 1), output: N x 1 
    Quadratic = @(x,beta) QuadraticForm(x,beta{1}) + sum(bsxfun(@times,x,beta{2}),2) + beta{3};
    
    function SMC = RunSMC(Twisting) 
       % Controlled SMC with a given policy
           % Input arguments: Twisting (struct)
           %                  (fields)
           %                    A (T+1 x d x d)
           %                    b (T+1 x d)
           %                    c (T+1 x 1)
           %                    G (T x d x d)
           %                    LogDetG (T x 1)
           % Output arguments: SMC (struct)
           %                   (fields)
           %                    Trajectory (N x d x T+1)
           %                    Ancestry (N x T+1)
           %                    LogIncrementalWeight (N x T+1)
           %                    ESS (1 x T+1)
           %                    LogNormConst (1 x T+1)
           %                    Trajectory (N x d x T+1)
           %                    ShiftTransitionMean (N x d x T)
           %                    QuadForm (N x T)
           %

        % Pre-allocate
        SMC = struct();
        SMC.Trajectory = zeros(N,d,T+1);
        SMC.Ancestry = zeros(N,T+1);
        SMC.LogIncrementalWeight = zeros(N,T+1);
        SMC.NormalisedWeights = zeros(N,T+1);
        SMC.ESS = zeros(1,T+1);
        SMC.LogNormConst = zeros(1,T+1);  
        SMC.ShiftTransitionMean = zeros(N,d,T);
        SMC.QuadForm = zeros(N,T);
        
        % Initialise
        TwistedInitialPrecision = Initial.Precision + 2 * squeeze(Twisting.A(1,:,:));
        TwistedInitialCov = inv(TwistedInitialPrecision);
        TwistedInitialMean = ( Initial.MeanPrecision - Twisting.b(1,:) ) * TwistedInitialCov;
        X = mvnrnd(TwistedInitialMean, TwistedInitialCov, N);
        SMC.Trajectory(:,:,1) = X;
        
        % Compute first transition kernel and log normalising function
        TransitionMean = Transition.Mean(1,X);
        ShiftTransitionMean = bsxfun(@minus, TransitionMean * Transition.Ginv, Twisting.b(2,:));
        SMC.ShiftTransitionMean(:,:,1) = ShiftTransitionMean; 
        TwistedG = squeeze(Twisting.G(1,:,:));
        TwistedTransitionMean = ShiftTransitionMean * TwistedG;
        QuadForm = 0.5 * QuadraticForm(ShiftTransitionMean, TwistedG); 
        SMC.QuadForm(:,1) = QuadForm;
        LogNormalisingFunc = 0.5 * Twisting.LogDetG(1) - 0.5 * Transition.LogDetG ... 
                           - 0.5 * QuadraticForm(TransitionMean, Transition.Ginv) - Twisting.c(2) ...
                           + QuadForm;        
        
        % Compute first weight function        
        TwistedInitialBeta = {squeeze(Twisting.A(1,:,:)), Twisting.b(1,:), Twisting.c(1)};
        LogNormalisingConst = 0.5 * log(det(TwistedInitialCov)) - 0.5 * Initial.LogDetCov ...
                            - 0.5 * Initial.QuadForm - Twisting.c(1) ... 
                            + 0.5 * sum((TwistedInitialMean * TwistedInitialPrecision) .* TwistedInitialMean);

        LogWeights = Observation.LogDensity(1,X) + LogNormalisingConst ... 
                   + LogNormalisingFunc  + Quadratic(X,TwistedInitialBeta);
        
        % Normalize weights, compute ESS and normalizing constant
        MaxLogWeight = max(LogWeights);
        Weights = exp(LogWeights - MaxLogWeight);
        SMC.LogIncrementalWeight(:,1) = LogWeights - LogNormalisingConst;
        NormalisedWeights = Weights / sum(Weights);
        SMC.NormalisedWeights(:,1) = NormalisedWeights;
        SMC.ESS(1) = 1 / sum(NormalisedWeights.^2);
        LogRatioNormConst = log(mean(Weights)) + MaxLogWeight;
        SMC.LogNormConst(1) = LogRatioNormConst;
                
        % Resample
        Ancestor = Systematic_Resampling(NormalisedWeights);
        SMC.Ancestry(:,1) = Ancestor;
        TwistedTransitionMean = TwistedTransitionMean(Ancestor,:);

        
        for t = 1:T
            % Move 
            X = mvnrnd(TwistedTransitionMean, TwistedG);
            SMC.Trajectory(:,:,t+1) = X;
            
            if (t < T)
                % Compute next transition kernel and log normalising function
                TransitionMean = Transition.Mean(t+1,X);
                ShiftTransitionMean = bsxfun(@minus, TransitionMean * Transition.Ginv, Twisting.b(t+2,:));
                SMC.ShiftTransitionMean(:,:,t+1) = ShiftTransitionMean; 
                TwistedG = squeeze(Twisting.G(t+1,:,:));
                TwistedTransitionMean = ShiftTransitionMean * TwistedG;
                QuadForm = 0.5 * QuadraticForm(ShiftTransitionMean, TwistedG);
                SMC.QuadForm(:,t+1) = QuadForm;
                LogNormalisingFunc = 0.5 * Twisting.LogDetG(t+1) - 0.5 * Transition.LogDetG ... 
                               - 0.5 * QuadraticForm(TransitionMean, Transition.Ginv) - Twisting.c(t+2) ...
                               + QuadForm;
            
            else
                LogNormalisingFunc = zeros(N,1);
            end
            
            % Compute incremental weights
            TwistedTransitionBeta = {squeeze(Twisting.A(t+1,:,:)), Twisting.b(t+1,:), Twisting.c(t+1)};
            LogWeights = Observation.LogDensity(t+1,X) + LogNormalisingFunc + Quadratic(X,TwistedTransitionBeta);
            SMC.LogIncrementalWeight(:,t+1) = LogWeights;
            
            % Normalize weights, compute ESS and normalizing constant
            MaxLogWeight = max(LogWeights);
            Weights = exp(LogWeights - MaxLogWeight);
            NormalisedWeights = Weights / sum(Weights);
            SMC.NormalisedWeights(:,t+1) = NormalisedWeights;
            SMC.ESS(t+1) = 1 / sum(NormalisedWeights.^2);
            LogRatioNormConst = LogRatioNormConst + log(mean(Weights)) + MaxLogWeight;
            SMC.LogNormConst(t+1) = LogRatioNormConst;   
            
            % Resampling
            Ancestor = Systematic_Resampling(NormalisedWeights);
            SMC.Ancestry(:,t+1) = Ancestor;
            TwistedTransitionMean = TwistedTransitionMean(Ancestor,:);

        end
          
    end
    
   % Initialise with bootstrap particle filter
    CurrentPolicy = struct();
    RepG = zeros(1,d,d);
    RepG(1,:,:) = Transition.G;
    CurrentPolicy.G = repmat(RepG,[T 1 1]);
    CurrentPolicy.LogDetG = Transition.LogDetG * ones(T,1);
    CurrentPolicy.A = zeros(T+1,d,d);
    CurrentPolicy.b = zeros(T+1,d);
    CurrentPolicy.c = zeros(T+1,1);    
    
    OutputSMC = RunSMC(CurrentPolicy); 
    cSMC = struct();
    cSMC.Trajectory = squeeze(OutputSMC.Trajectory);
    cSMC.SmoothingWeights = zeros(N,T+1);
    LogSmoothingWeights = log( OutputSMC.NormalisedWeights(:,end) );
    cSMC.SmoothingWeights(:,end) = OutputSMC.NormalisedWeights(:,end);
    
    for s = T:-1:1
        LogFilteringWeights = log( OutputSMC.NormalisedWeights(:,s) );     
        LogPredictiveMeasure = zeros(N,N);        
        for n = 1:N
            LogPredictiveMeasure(:,n) = LogFilteringWeights + Transition.LogDensity( OutputSMC.Trajectory(:,1,s) , OutputSMC.Trajectory(n,1,s+1) );
        end
        MaxLogMeasure = max(max(LogPredictiveMeasure));
        LogSumPredictive = log( sum( exp(LogPredictiveMeasure - MaxLogMeasure) , 1) );
        
        LogJointMeasure = zeros(N,N);
        for j = 1:N
            LogJointMeasure(:,j) = LogSmoothingWeights(j) + LogPredictiveMeasure(:,j) - MaxLogMeasure - LogSumPredictive(j);
        end
        MaxLogJointMeasure = max(max(LogJointMeasure));
        
        LogSmoothingWeights = log( sum( exp(LogJointMeasure - MaxLogJointMeasure) , 2) ) + MaxLogJointMeasure;
        cSMC.SmoothingWeights(:,s) = exp(LogSmoothingWeights); 
        
    end
    
end

