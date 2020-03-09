function [V] = Volume(X)
    [r,~] = size(unique(X,'rows'));
    if r < 3
        V = 0;
    else
        [~,V] = convhull(X(:,1),X(:,2));
    end
    
end

