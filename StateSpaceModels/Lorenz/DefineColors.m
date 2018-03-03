function colors = DefineColors(I)
% Define colors for each iteration
%   Input argument: I (1 x 1) 
%   Output argument: colors (cell-array)

    colour_teal = [18 150 155] ./ 255;
    colour_lightgreen = [94 250 81] ./ 255;
    colour_lightblue = [8 180 238] ./ 255;
    colour_darkblue = [1 17 181] ./ 255;
    colour_gold = [225 228 42] ./ 255;
    colour_peach = [251 111 66] ./ 255;
    colour_black = [0 0 0];
    colour_magenta = [1 0 1];
    colour_purple = [121 84 233] ./ 255;

    switch I 
        case 0 
            colors = {colour_darkblue};
        case 1 
            colors = {colour_darkblue, colour_black};
        case 2
            colors = {colour_darkblue, colour_black, colour_teal}; 
        case 3
            colors = {colour_darkblue, colour_black, colour_teal, colour_magenta};
        case 4
            colors = {colour_darkblue, colour_black, colour_teal, colour_magenta, colour_peach}; 
        case 5
            colors = {colour_darkblue, colour_black, colour_teal, colour_magenta, colour_peach, colour_gold};
        case 6
            colors = {colour_darkblue, colour_black, colour_teal, colour_magenta, colour_peach, ... 
                colour_gold, colour_lightblue};
        case 7
            colors = {colour_darkblue, colour_black, colour_teal, colour_magenta, colour_peach, ...
                colour_gold, colour_lightblue, colour_lightgreen}; 
        case 8
            colors = {colour_darkblue, colour_black, colour_teal, colour_magenta, colour_peach, ... 
                colour_gold, colour_lightblue, colour_lightgreen, colour_purple}; 
    end


end

