function V = PCoVol(X)
    % PCo Volumne in N dimensions. Rows of X must be PC scores,
    % columns of X must be variables (e.g:       | var 1 | var 2 | var 3 |
    %                                     taxa 1 | 1.293 | ...
    %                                     taxa 2 | 0.370 | ...
    %                                     taxa 3 | 2.100 | ...
    
    
    unq = unique(X,'rows');
    [L,~] = size(unq);
    if L < 2
        V = nan;
    else
        [~,scr,lat] = pca(X);
        sumL = sum(lat);
        tot = 0;
        for i = 1:length(lat)
            ex = lat(i)/sumL;
            if tot < 0.95
                tot = tot + ex;
                if tot > 0.95
                    N = i;
                end
            end
        end
        scr = scr(:,1:N);

        if N < 2
            V = max(scr(:,1)) - min(scr(:,1));
        elseif L < 3
            V = nan;
        elseif N < 4
            [~,V] = convhull(scr);
        else
            [~,V] = convhulln(scr);
        end
        V = nthroot(V,N);
    end
end

