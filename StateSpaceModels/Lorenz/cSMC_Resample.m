function cSMC = cSMC_Resample(Parameters, Initial, Transition, Observation, InitialPolicy)
% Controlled sequential Monte Carlo
    % Input arguments: Parameters (struct)
    %                  Initial (struct)
    %                  Transition (struct)
    %                  Observation (struct)
    %                  InitialPolicy (struct)
    % Output arguments: cSMC (cell-array of size I+1 x 2)

    N = Parameters.Particles;
    d = Parameters.Dimension;
    T = Parameters.Time;
    I = Parameters.Iterations;
    Approximator = Parameters.Approximator;
    
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
        SMC.ESS = zeros(1,T+1);
        SMC.LogNormConst = zeros(1,T+1);  
        SMC.ShiftTransitionMean = zeros(N,d,T);
        SMC.QuadForm = zeros(N,T);
        
        % Initialise
        TwistedInitialPrecision = Initial.Precision + 2 * squeeze(Twisting.A(1,:,:));
        TwistedInitialCov = inv(TwistedInitialPrecision);
        TwistedInitialMean = ( Initial.MeanPrecision - Twisting.b(1,:) ) * TwistedInitialCov;
%         TwistedInitialCov = EnsureSPD(TwistedInitialCov); % ensure positive definite
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
        SMC.ESS(1) = 1 / sum(NormalisedWeights.^2);
        LogRatioNormConst = log(mean(Weights)) + MaxLogWeight;
        SMC.LogNormConst(1) = LogRatioNormConst;
                
        % Resample
        Ancestor = Systematic_Resampling(NormalisedWeights);
        SMC.Ancestry(:,1) = Ancestor;
        TwistedTransitionMean = TwistedTransitionMean(Ancestor,:);

        
        for t = 1:T
            % Move 
%             TwistedG = EnsureSPD(TwistedG); % ensure positive definite
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
            SMC.ESS(t+1) = 1 / sum(NormalisedWeights.^2);
            LogRatioNormConst = LogRatioNormConst + log(mean(Weights)) + MaxLogWeight;
            SMC.LogNormConst(t+1) = LogRatioNormConst;   
            
            % Resampling
            Ancestor = Systematic_Resampling(NormalisedWeights);
            SMC.Ancestry(:,t+1) = Ancestor;
            TwistedTransitionMean = TwistedTransitionMean(Ancestor,:);

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
    Indexc = 1;
    
    function NewPolicy = ApproxDP(SMC,Policy)
        % Approximate dynamic programming 
            % Input argument: SMC (struct)
            %                 (fields)
            %                   LogIncrementalWeight (N x T+1)     
            %                   Trajectory (N x d x T+1)
            %                   ShiftTransitionMean (N x d x T)
            %                   QuadForm (N x T)
            %                 Policy (struct)
            %                 (fields)
            %                   A (T+1 x d x d)
            %                   b (T+1 x d)
            %                   c (T+1 x 1)
            % Output argument: NewPolicy (struct)
            %                  (fields)
            %                   A (T+1 x d x d)
            %                   b (T+1 x d)
            %                   c (T+1 x 1)
            %                   G (T x d x d)
            %                   LogDetG (T x 1)
            
        % Pre-allocate
        NewPolicy = struct();
        NewPolicy.A = zeros(T+1,d,d);
        NewPolicy.b = zeros(T+1,d);
        NewPolicy.c = zeros(T+1,1);
        NewPolicy.G = zeros(T,d,d);
        NewPolicy.LogDetG = zeros(T,1);
        
        % Initialise
        LogNormalisingFunc = zeros(N,1);
        
        for t = (T+1):-1:1
            % Compute value at particle locations
            Value = - SMC.LogIncrementalWeight(:,t) - LogNormalisingFunc;
            
            % Linear least squares
            DesignMatrix = x2fx(squeeze(SMC.Trajectory(:,:,t)),Approximator);
            BetaHat = (DesignMatrix' * DesignMatrix) \ (DesignMatrix' * Value); % (d^2/2 + 3d/2 + 1) x 1
            
            % Store parameters
            if ( d == 1 || strcmp(Approximator,'purequadratic') )
                AHatVec = BetaHat(IndexAA);
                if (t == 1)
                    PDconstraint = - 0.5 * diag(Initial.Precision) - diag(squeeze(Policy.A(1,:,:)));    
                else
                    PDconstraint = - 0.5 * diag(Transition.Ginv) - diag(squeeze(Policy.A(t,:,:)));    
                end
                IndexConstraint = (AHatVec < PDconstraint);
                if ( any(IndexConstraint) )
                    disp(['Projecting at time ' num2str(t)])
                end
                AHatVec(IndexConstraint) = PDconstraint(IndexConstraint) + 0.001;
                AHat = diag(AHatVec);
%                 AHat = diag(BetaHat(IndexAA));
            else
                AHat = diag(BetaHat(IndexAA)) + 0.5 * squareform(BetaHat(IndexA));
            end
            bHat = BetaHat(Indexb)';
            cHat = BetaHat(Indexc);
            
            UpdatedA = squeeze(Policy.A(t,:,:)) + AHat;
            
            NewPolicy.A(t,:,:) = UpdatedA;
            NewPolicy.b(t,:) = Policy.b(t,:) + bHat;
            NewPolicy.c(t) = Policy.c(t) + cHat; 
            
            if (t > 1)
                G = inv( Transition.Ginv + 2 * UpdatedA );
                LogDetG = log(det(G));
                NewPolicy.G(t-1,:,:) = G;
                NewPolicy.LogDetG(t-1) = LogDetG;

                % Evaluate next log normalising function
                NewShiftTransitionMean = bsxfun(@minus, SMC.ShiftTransitionMean(:,:,t-1) , bHat);
                LogNormalisingFunc = 0.5 * LogDetG - 0.5 * Policy.LogDetG(t-1) - cHat ... 
                                   + 0.5 * QuadraticForm(NewShiftTransitionMean, G) - SMC.QuadForm(:,t-1);         
            end
            
        end
              
    end

    % Pre-allocate
    cSMC = cell(I+1,2); % store value functions (left) and SMC output (right)
    CurrentPolicy = InitialPolicy;
    cSMC{1,1} = CurrentPolicy;
    OutputSMC = RunSMC(CurrentPolicy); 
    cSMC{1,2} = OutputSMC;
    
    switch Parameters.Adaptive
        case 'yes'
            for i = 1:I
                UpdatedPolicy = ApproxDP(OutputSMC,CurrentPolicy);         
                cSMC{i+1,1} = UpdatedPolicy;
                CurrentPolicy = UpdatedPolicy;

                OutputSMC = RunSMC(CurrentPolicy); 
                cSMC{i+1,2} = OutputSMC;
                if and( min(OutputSMC.ESS) * 100 / N > 90, i<I)
                    cSMC(i+2:I+1,:) = [];
                    break
                end
            end
        case 'no'
            % Iterate 
            for i = 1:I
                UpdatedPolicy = ApproxDP(OutputSMC,CurrentPolicy);         
                cSMC{i+1,1} = UpdatedPolicy;
                CurrentPolicy = UpdatedPolicy;

                OutputSMC = RunSMC(CurrentPolicy); 
                cSMC{i+1,2} = OutputSMC;               
            end
    end

end

