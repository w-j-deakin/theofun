function V = Vol3D(X)
    
    
    [r,c] = size(X);
    
    if c > 3
        x = X(:,1:3);
    elseif c < 3
        x = zeros(r,3);
        for i = 1:c
            x(:,c) = X(:,c);
        end
    else
        x = X;
    end
    
    unq = unique(x,'rows');
    [L,~] = size(unq);
    if L == 1
        V = 0;
    elseif L == 2
        v = unq(1,:) - unq(2,:);
        V = norm(v);
    elseif L == 3
        v = x(2,:) - x(1,:);
        w = x(3,:) - x(1,:);
        V = sqrt(0.5*norm(cross(v,w)));
    else
        [~,v] = convhull(x(:,1),x(:,2),x(:,3));
        V = nthroot(v,3);
    end
    
end

