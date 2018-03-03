function Output = EvaluatePolicy(t,grid,Policy)
% Evaluate twisting function 

    SummandPolicy = bsxfun(@plus, log(Policy.Weights(:,t))', - Policy.Bandwidth(t) * bsxfun(@minus, grid, Policy.Knots(:,t)').^2); % N x M
    MaxSummandPolicy = max(max(SummandPolicy));
    LogPolicy = log( sum( exp(SummandPolicy - MaxSummandPolicy), 2) ) + MaxSummandPolicy + Policy.MaxLogPolicyValue(t);
    Output = exp(LogPolicy);

end

