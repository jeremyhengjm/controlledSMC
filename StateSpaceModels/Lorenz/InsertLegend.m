function h = InsertLegend(I)
% Include legend for each iteration 
%   Input argument: I (1 x 1)
%   Output argument: h (function handle)

    switch I
        case 0 
            h = legend('Uncontrolled');
        case 1
             h = legend('Uncontrolled','Iteration 1');
        case 2
            h = legend('Uncontrolled','Iteration 1','Iteration 2','Location','Best');   
        case 3 
            h = legend('Uncontrolled','Iteration 1','Iteration 2','Iteration 3','Location','Best');   
        case 4
            h = legend('Uncontrolled','Iteration 1','Iteration 2','Iteration 3','Iteration 4','Location','Best');     
        case 5
            h = legend('Uncontrolled','Iteration 1','Iteration 2','Iteration 3','Iteration 4', ... 
                        'Iteration 5');
        case 6
            h = legend('Uncontrolled','Iteration 1','Iteration 2','Iteration 3','Iteration 4', ... 
                        'Iteration 5','Iteration 6','Location','Best'); 
        case 7
            h = legend('Uncontrolled','Iteration 1','Iteration 2','Iteration 3','Iteration 4', ... 
                        'Iteration 5','Iteration 6','Iteration 7','Location','Best'); 
        case 8
            h = legend('Uncontrolled','Iteration 1','Iteration 2','Iteration 3','Iteration 4', ... 
                        'Iteration 5','Iteration 6','Iteration 7','Iteration 8','Location','Best'); 
    end


end

