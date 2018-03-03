function Ancestor = Systematic_Resampling(NormalisedWeights)
% Systematic resampling
    % Input argument: NormalisedWeights (N x 1)
    % Output argument: Ancestor (N x 1)
    
    N = length(NormalisedWeights);
    Ancestor = zeros(N,1);
    NormalisedWeights = N * NormalisedWeights;
    j = 1;
    csw = NormalisedWeights(1);
    U = rand(1);
    for n = 1:N
        while (csw < U)
            j = j + 1;
            csw = csw + NormalisedWeights(j);
        end
        Ancestor(n) = j;
        U = U + 1;
    end

end

