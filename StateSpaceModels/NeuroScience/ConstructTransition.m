function Transition = ConstructTransition(theta)
% Construct transition kernel f

    Transition = struct();
    Transition.Mean = @(t,x) x * theta(1);
    Transition.G = theta(2);
    Transition.Ginv = 1 / theta(2);
    Transition.LogDetG = log(theta(2));

end

