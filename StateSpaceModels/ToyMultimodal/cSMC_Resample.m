function cSMC = cSMC_Resample( Parameters, Initial, Transition, Observation )
% Controlled Sequential Monte Carlo
    % Input arguments: Parameters (struct)
    %                  Initial (struct)
    %                  Transition (struct)
    %                  Observation (struct)
    % Output arguments: cSMC (cell-array of size I+1 x 2)

    N = Parameters.Particles;
    M = Parameters.Knots;
    T = Parameters.Time;
    I = Parameters.Iterations;
    BandwidthFactor = Parameters.BandwidthFactor;

    RBF = @(beta,distancesq) exp(- beta * distancesq);
    
    function SMC = RunSMC(Policy)
    % Controlled SMC with a given policy 
    % Input argument: Policy (struct)
    %                 (fields)
    %                   Bandwidth (1 x T+1)
    %                   Weights (M x T+1)
    %                   Knots (M x T+1)
    % Output argument: SMC (struct()
    %                   (fields)
    %                    Trajectory (N x T+1)
    %                    Ancestry (N x T+1)
    %                    LogIncrementalWeight (N x T+1)
    %                    ESS (1 x T+1)
    %                    LogNormConst (1 x T+1)
    %                    LogObservationDensity (N x T+1)
    %                    TransitionMean (N x T)
    
        % Pre-allocate
        SMC = struct();
        SMC.Trajectory = zeros(N,T+1);
        SMC.Ancestry = zeros(N,T+1);
        SMC.LogIncrementalWeight = zeros(N,T+1);
        SMC.ESS = zeros(1,T+1);
        SMC.LogNormConst = zeros(1,T+1);     
        SMC.LogObservationDensity = zeros(N,T+1);
        SMC.TransitionMean = zeros(N,T);

        % Initialise
        TwistedInitialPrecision = Initial.Precision + 2 * Policy.Bandwidth(1);
        TwistedInitialVar = 1 / TwistedInitialPrecision; % 1 x 1
        TwistedInitialMean = TwistedInitialVar * ( Initial.MeanPrecision + ... 
                                2 * Policy.Bandwidth(1) * Policy.Knots(:,1) ); % M x 1
                            
        LogPolicyWeights = log(Policy.Weights(:,1)) - Policy.Bandwidth(1) * Policy.Knots(:,1).^2 + ... 
            0.5 * TwistedInitialPrecision * TwistedInitialMean.^2; % M x 1
        MaxLogPolicyWeights = max(LogPolicyWeights);
        PolicyWeights = exp(LogPolicyWeights - MaxLogPolicyWeights);
        SumPolicyWeights = sum(PolicyWeights);
        NormalisedPolicyWeights = PolicyWeights / SumPolicyWeights;  
        SelectComponents = DiscreteSample(NormalisedPolicyWeights, N); % select mixture components
        X = TwistedInitialMean(SelectComponents) + sqrt(TwistedInitialVar) * randn(N,1);
        SMC.Trajectory(:,1) = X;
        
        % Compute log normalising constant
        LogSum = log(SumPolicyWeights) + MaxLogPolicyWeights; % 1 x 1 
        LogNormalisingConst = 0.5 * log(TwistedInitialVar) - 0.5 * log(Initial.Var) ...
                            - 0.5 * Initial.Mean^2 * Initial.Precision + LogSum + Policy.MaxLogPolicyValue(1); % 1 x 1
                        
        % Compute first log normalising function
        TwistedTransitionPrecision = Transition.Ginv + 2 * Policy.Bandwidth(2);
        TwistedTransitionVar = 1 / TwistedTransitionPrecision;
        TransitionMean = Transition.Mean(1,X); % N x 1
        SMC.TransitionMean(:,1) = TransitionMean;
        TwistedTransitionMean = TwistedTransitionVar * bsxfun(@plus, Transition.Ginv * TransitionMean, ... 
                                    2 * Policy.Bandwidth(2) * Policy.Knots(:,2)' ); % N x M
                                
        LogPolicyWeights = bsxfun(@plus, 0.5 * TwistedTransitionPrecision * TwistedTransitionMean.^2, ...
            log(Policy.Weights(:,2))' - Policy.Bandwidth(2) * Policy.Knots(:,2).^2' ); % N x M
        MaxLogPolicyWeights = max(max(LogPolicyWeights)); 
        PolicyWeights = exp(LogPolicyWeights - MaxLogPolicyWeights); % N x M
        SumPolicyWeights = sum(PolicyWeights, 2); % N x 1

        LogSumTwisted = log(SumPolicyWeights) + MaxLogPolicyWeights; % N x 1
        
        LogNormalisingFunc = 0.5 * log(TwistedTransitionVar) - 0.5 * Transition.LogG  ... 
                           - 0.5 * Transition.Ginv * TransitionMean.^2 + LogSumTwisted + Policy.MaxLogPolicyValue(2); % N x 1
       
        % Compute first weight function 
        LogObservationDensity = Observation.LogDensity(1,X); % N x 1
        SMC.LogObservationDensity(:,1) = LogObservationDensity;
        
        SummandPolicy = bsxfun(@plus, log(Policy.Weights(:,1))', - Policy.Bandwidth(1) * bsxfun(@minus, X, Policy.Knots(:,1)').^2); % N x M
        MaxSummandPolicy = max(max(SummandPolicy));
        LogPolicy = log( sum( exp(SummandPolicy - MaxSummandPolicy), 2) ) + MaxSummandPolicy + Policy.MaxLogPolicyValue(1);
        
        LogWeights = LogObservationDensity + LogNormalisingConst ... 
                   + LogNormalisingFunc  - LogPolicy;
               
        % Normalize weights, compute ESS and normalizing constant
        MaxLogWeight = max(LogWeights);
        Weights = exp(LogWeights - MaxLogWeight);
        NormalisedWeights = Weights / sum(Weights);
        SMC.ESS(1) = 1 / sum(NormalisedWeights.^2);
        LogRatioNormConst = log(mean(Weights)) + MaxLogWeight;
        SMC.LogNormConst(1) = LogRatioNormConst;
        SMC.LogIncrementalWeight(:,1) = LogWeights - LogNormalisingConst;
                
        % Resample
        Ancestor = Systematic_Resampling(NormalisedWeights);
        SMC.Ancestry(:,1) = Ancestor;
            
        % Iterate
        for t = 1:T
            NormalisedPolicyWeights = bsxfun(@rdivide, PolicyWeights, SumPolicyWeights); % N x M        
            
            SelectComponents = zeros(1,N);
            for n = 1:N
                SelectComponents(n) = DiscreteSample(NormalisedPolicyWeights(n,:), 1);
            end
            
            IndexAncestorSelect = sub2ind([N M], Ancestor', SelectComponents);
            X = TwistedTransitionMean(IndexAncestorSelect)' + sqrt(TwistedTransitionVar) * randn(N,1);
            SMC.Trajectory(:,t+1) = X;
            
            if (t < T)
                % Compute next log normalising function
                TwistedTransitionPrecision = Transition.Ginv + 2 * Policy.Bandwidth(t+2);
                TwistedTransitionVar = 1 / TwistedTransitionPrecision;
                TransitionMean = Transition.Mean(t+1,X); % N x 1
                SMC.TransitionMean(:,t+1) = TransitionMean;
                TwistedTransitionMean = TwistedTransitionVar * bsxfun(@plus, Transition.Ginv * TransitionMean, ... 
                                        2 * Policy.Bandwidth(t+2) * Policy.Knots(:,t+2)' ); % N x M

                LogPolicyWeights = bsxfun(@plus, 0.5 * TwistedTransitionPrecision * TwistedTransitionMean.^2, ...
                        log(Policy.Weights(:,t+2))' - Policy.Bandwidth(t+2) * Policy.Knots(:,t+2).^2' ); % N x M
                MaxLogPolicyWeights = max(max(LogPolicyWeights)); 
                PolicyWeights = exp(LogPolicyWeights - MaxLogPolicyWeights); % N x M
                SumPolicyWeights = sum(PolicyWeights, 2); % N x 1
        
                LogSumTwisted = log(SumPolicyWeights) + MaxLogPolicyWeights; % N x 1
        
                LogNormalisingFunc = 0.5 * log(TwistedTransitionVar) - 0.5 * Transition.LogG  ... 
                                   - 0.5 * Transition.Ginv * TransitionMean.^2 + LogSumTwisted + Policy.MaxLogPolicyValue(t+2); % N x 1                 
            else
                LogNormalisingFunc = zeros(N,1);
            end
            
            % Compute incremental weights
            LogObservationDensity = Observation.LogDensity(t+1,X); % N x 1
            SMC.LogObservationDensity(:,t+1) = LogObservationDensity;
        
            SummandPolicy = bsxfun(@plus, log(Policy.Weights(:,t+1))', - Policy.Bandwidth(t+1) * bsxfun(@minus, X, Policy.Knots(:,t+1)').^2); % N x M
            MaxSummandPolicy = max(max(SummandPolicy));
            LogPolicy = log( sum( exp(SummandPolicy - MaxSummandPolicy), 2) ) + MaxSummandPolicy + Policy.MaxLogPolicyValue(t+1);
        
            LogWeights = LogObservationDensity + LogNormalisingFunc - LogPolicy;
            SMC.LogIncrementalWeight(:,t+1) = LogWeights;    
            
            % Normalize weights, compute ESS and normalizing constant
            MaxLogWeight = max(LogWeights);
            Weights = exp(LogWeights - MaxLogWeight);
            NormalisedWeights = Weights / sum(Weights);
            SMC.ESS(t+1) = 1 / sum(NormalisedWeights.^2);
            LogRatioNormConst = LogRatioNormConst + log(mean(Weights)) + MaxLogWeight;
            SMC.LogNormConst(t+1) = LogRatioNormConst;   
            
            % Resampling
            Ancestor = Systematic_Resampling(NormalisedWeights);
            SMC.Ancestry(:,t+1) = Ancestor;

        end    
        
    end

    function NewPolicy = ApproxDP(SMC)
    % Approximate dynamic programming with radial basis functions
        % Input argument: SMC (struct)
        %                 (fields)
        %                   LogObservationDensity (N x T+1)
        %                   Trajectory (N x T+1)
        %                   TransitionMean (N x T)
        % Output argument: NewPolicy (struct)
        %                  (fields)
        %                   Bandwidth (1 x T+1)
        %                   Weights (M x T+1)
        %                   Knots (M x T+1)

        % Pre-allocate
        NewPolicy = struct();
        NewPolicy.Bandwidth = zeros(1,T+1);
        NewPolicy.Weights = zeros(M,T+1);
        NewPolicy.Knots = zeros(M,T+1);
        NewPolicy.MaxLogPolicyValue = zeros(1,T+1);
        
        % Initialise
        LogNormalisingFunc = zeros(N,1);
        
        switch Parameters.KnotsType
            case 'particles'
                for t = (T+1):-1:1
                    % Compute value at particle locations
                    LogPolicyValue = SMC.LogObservationDensity(:,t) + LogNormalisingFunc;
                    MaxLogPolicyValue = max(LogPolicyValue);
                    NewPolicy.MaxLogPolicyValue(t) = MaxLogPolicyValue;
                    PolicyValue = exp(LogPolicyValue - MaxLogPolicyValue);

                    % Fit radial basis function
                    Particles = SMC.Trajectory(:,t);
                    Bandwidth = std(Particles) * BandwidthFactor; 
                    NewPolicy.Bandwidth(t) = Bandwidth; 
                    PairwiseDistancesq = bsxfun(@minus, Particles, Particles').^2;
                    FittedWeights = lsqnonneg( RBF(Bandwidth, PairwiseDistancesq), PolicyValue );

                    % Remove components with low weights
                    [SortedWeights, IndexWeights] = sort(FittedWeights, 'descend');
                    PolicyWeights = SortedWeights(1:M);
                    NewPolicy.Weights(:,t) = PolicyWeights;
                    Knots = Particles(IndexWeights(1:M)); % M x 1
                    NewPolicy.Knots(:,t) = Knots;

                    % Compute next log normalising function
                    if (t > 1)
                        TwistedTransitionPrecision = Transition.Ginv + 2 * Bandwidth;
                        TwistedTransitionVar = 1 / TwistedTransitionPrecision;            
                        TransitionMean = SMC.TransitionMean(:,t-1); % N x 1

                        TwistedTransitionMean = TwistedTransitionVar * bsxfun(@plus, Transition.Ginv * TransitionMean, ... 
                                                2 * Bandwidth * Knots' ); % N x M

                        SummandTwisted = bsxfun(@plus, 0.5 * TwistedTransitionPrecision * TwistedTransitionMean.^2, ...
                                log(PolicyWeights)' - Bandwidth * Knots.^2' ); % N x M                    

                        MaxSummandTwisted = max(max(SummandTwisted)); 
                        LogSumTwisted = log( sum( exp(SummandTwisted - MaxSummandTwisted), 2) ) + MaxSummandTwisted; % N x 1

                        LogNormalisingFunc = 0.5 * log(TwistedTransitionVar) - 0.5 * Transition.LogG  ... 
                                           - 0.5 * Transition.Ginv * TransitionMean.^2 + LogSumTwisted + MaxLogPolicyValue; % N x 1             

                    end

                end
                
            case 'grid'
                Grid = Parameters.Grid;
                for t = (T+1):-1:1
                    % Compute value at particle locations
                    LogPolicyValue = Observation.LogDensity(t,Grid) + LogNormalisingFunc;
                    MaxLogPolicyValue = max(LogPolicyValue);
                    NewPolicy.MaxLogPolicyValue(t) = MaxLogPolicyValue;
                    PolicyValue = exp(LogPolicyValue - MaxLogPolicyValue);

                    % Fit radial basis function                    
                    Bandwidth = (Grid(end) - Grid(1)) * BandwidthFactor; 
                    NewPolicy.Bandwidth(t) = Bandwidth; 
                    PairwiseDistancesq = bsxfun(@minus, Grid, Grid').^2;
                    FittedWeights = lsqnonneg( RBF(Bandwidth, PairwiseDistancesq), PolicyValue );

                    % Remove components with low weights
                    [SortedWeights, IndexWeights] = sort(FittedWeights, 'descend');
                    PolicyWeights = SortedWeights(1:M);
                    NewPolicy.Weights(:,t) = PolicyWeights;
                    Knots = Grid(IndexWeights(1:M)); % M x 1
                    NewPolicy.Knots(:,t) = Knots;

                    % Compute next log normalising function
                    if (t > 1)
                        TwistedTransitionPrecision = Transition.Ginv + 2 * Bandwidth;
                        TwistedTransitionVar = 1 / TwistedTransitionPrecision; 
                        TransitionMean = Transition.Mean(t-1,Grid);                        

                        TwistedTransitionMean = TwistedTransitionVar * bsxfun(@plus, Transition.Ginv * TransitionMean, ... 
                                                2 * Bandwidth * Knots' ); % N x M

                        SummandTwisted = bsxfun(@plus, 0.5 * TwistedTransitionPrecision * TwistedTransitionMean.^2, ...
                                log(PolicyWeights)' - Bandwidth * Knots.^2' ); % N x M                    

                        MaxSummandTwisted = max(max(SummandTwisted)); 
                        LogSumTwisted = log( sum( exp(SummandTwisted - MaxSummandTwisted), 2) ) + MaxSummandTwisted; % N x 1

                        LogNormalisingFunc = 0.5 * log(TwistedTransitionVar) - 0.5 * Transition.LogG  ... 
                                           - 0.5 * Transition.Ginv * TransitionMean.^2 + LogSumTwisted + MaxLogPolicyValue; % N x 1             

                    end

                end
                
                
        end
        

    end

    % Pre-allocate
    cSMC = cell(I+1,2); % store value functions (left) and SMC output (right)
    
    % Initialise with bootstrap particle filter
    CurrentPolicy = struct();
    CurrentPolicy.Bandwidth = zeros(1,T+1);
    CurrentPolicy.Weights = ones(M,T+1) / M ;
    CurrentPolicy.Knots = zeros(M,T+1);
    CurrentPolicy.MaxLogPolicyValue = zeros(1,T+1);
    cSMC{1,1} = CurrentPolicy;
    OutputSMC = RunSMC(CurrentPolicy); 
    cSMC{1,2} = OutputSMC;
    
    % Iterate 
    for i = 1:I
        CurrentPolicy = ApproxDP(OutputSMC); 
        cSMC{i+1,1} = CurrentPolicy;
        
        OutputSMC = RunSMC(CurrentPolicy); 
        cSMC{i+1,2} = OutputSMC;
    end

end

