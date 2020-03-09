function sov = soVar(X)
    % Sum of Variances. Rows of X must be observations, columns of X must
    % be variables (e.g:       | var 1 | var 2 | var 3 |
    %                   taxa 1 | 1.293 | ...
    %                   taxa 2 | 0.370 | ...
    %                   taxa 3 | 2.100 | ...
    
    
    V = var(X,0,1);
    sov = sum(V);
    
end

