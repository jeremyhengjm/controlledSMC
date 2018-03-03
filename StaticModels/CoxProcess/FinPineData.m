% Process Finland pine tree data
Grid = 30;
G = linspace(0,1,Grid+1);

Counts = zeros(Grid,Grid);

for i = 1:Grid
    for j = 1:Grid
        logi_x = and(Like.x > G(i),Like.x < G(i+1));
        logi_y = and(Like.y > G(j),Like.y < G(j+1));
        Counts(i,j) = sum(logi_x .* logi_y);
    end
end

Counts = Counts';
Like.Obs = Counts(:); % convert matrix into vector (column-wise)

