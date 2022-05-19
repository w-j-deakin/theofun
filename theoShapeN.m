classdef theoShapeN
    
    % Individual theoretical shape object.
    
    properties
        v
        harm
        nH
        isIntersect
    end
    
    methods
        function obj = theoShapeN(score,H)
            [r,c] = size(score);
            if r > 1 && c > 1
                error('input score must be a vector');
            elseif c > 1
                score = score';
            end
            
            obj.v = score;
            
            [r,c] = size(H);
            if c == 4
                harm = H;
                nHarm = r;
            elseif r < 1 && c < 1
                error('invalid harmnic input');
            else
                nH = length(H);

                nHarm = (nH + 3)/4;

                harm = zeros(nHarm,4);
                harm(1,1) = 1;
                harm(1,4) = H(1);

                for i = 2:nH
                    idx1 = ceil(i/4) + 1;
                    idx2 = mod((i-1),4);
                    if idx2 == 0
                        idx2 = 4;
                        idx1 = idx1-1;
                    end
                    harm(idx1,idx2) = H(i);
                end
            end
            obj.nH = nHarm;
            obj.harm = harm;
            obj.isIntersect = selfIntersect(harm);
        end
        
        function tf = checkTSN(obj)
            tf = obj.isIntersect;
        end
        
        function [x,y] = draw(obj,n,inX,inY,s)
            
            H = obj.harm;
            
            [nHarm,~] = size(H);
            
            X = zeros(n,nHarm);
            Y = zeros(n,nHarm);
            
            x = zeros(n,1);
            y = zeros(n,1);
            
            for i = 1:n
                for j = 1:nHarm
                    X(i,j) = (H(j,1) * cos(j * (i*2*pi)/n)) + (H(j,2) * sin(j * (i*2*pi)/n));
                    Y(i,j) = (H(j,3) * cos(j * (i*2*pi)/n)) + (H(j,4) * sin(j * (i*2*pi)/n));
                end 
                x(i) = sum(X(i,:));
                y(i) = sum(Y(i,:));
            end
            
            if nargin > 2
                x = (x*s) + inX;
                y = (y*s) + inY;
            end
            
        end

        %% functional calculations

        % wing aspect ratio
        function AR = AspectRatio(obj,n)
            if (obj.isIntersect)
                AR = nan;
            else
                [x,y] = obj.draw(n,0,0,1);
                [~,baseID] = min(x);
                [~,tipID] = max(x);

                L = norm([(x(tipID) - x(baseID)) (y(tipID) - y(baseID))]);

                A = polyarea(x,y);
                AR = L*L/A;
            end
        end

        % moment of area

        function [i1x, i1y, i2x, i2y, i2z] = MoA(obj, n)
            if (obj.isIntersect)
                i1x = nan;
                i1y = nan;
                i2x = nan;
                i2y = nan;
                i2z = nan;
            else
                [x,y] = obj.draw(n,0,0,1);
                [~,baseID] = min(x);
                [~,tipID] = max(x);
                A = findAngle2D([1 0], [(x(tipID) - x(baseID)) (y(tipID) - y(baseID))]);
                x = x - x(baseID);
                y = y - y(baseID);

                points = Rotate([x y], A);

                A = polyarea(points(:,1),points(:,2));
                points = points / sqrt(A);

                preI = n;
                i1x = 0;
                i1y = 0;
                i2x = 0;
                i2y = 0;
                for i = 1:n
                    x0 = points(preI,1);
                    x1 = points(i,1);
                    y0 = points(preI,2);
                    y1 = points(i,2);
                    A = x0*y1 - x1*y0;
                    i1x = i1x + A * (y0 + y1);
                    i1y = i1y + A * (x0 + x1);
                    i2x = i2x + A * (y0*y0 + y0*y1 + y1*y1);
                    i2y = i2y + A * (x0*x0 + x0*x1 + x1*x1);
                    preI = i;
                end
                
                i1x = abs(i1x / 6);
                i1y = abs(i1y / 6);
                i2x = abs(i2x / 12);
                i2y = abs(i2y / 12);
                i2z = i2x + i2y;
            end
        end
        
    end
end

