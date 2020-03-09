function d = meanDist(X)
    % mean pairwise distance. Rows of X must be observations, columns of X must
    % be variables (e.g:       | var 1 | var 2 | var 3 |
    %                   taxa 1 | 1.293 | ...
    %                   taxa 2 | 0.370 | ...
    %                   taxa 3 | 2.100 | ...
    
    D = pdist(X);
    d = mean(D);
end

