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
                L = max(x) - min(x);
                disp(L);
                A = polyarea(x,y);
                disp(A);
                AR = L*L/A;
            end
        end

        % second moment of area
        function [ix, iy] = SMOA(obj, n)
            if (obj.isIntersect)
                iy = nan;
                ix = nan;
            else
                [x,y] = obj.draw(n,0,0,1);
                A = polyarea(x,y);
                x = x / sqrt(A);
                y = y / sqrt(A);
                preI = n;
                iy = 0;
                ix = 0;
                for i = 1:n
                    x0 = x(preI);
                    x1 = x(i);
                    y0 = y(preI);
                    y1 = y(i);
                    ix = ix + (((x0*y1)-(x1*y0))*((y0*y0)+(y0*y1)+(y1*y1)));
                    iy = iy + (((x0*y1)-(x1*y0))*((x0*x0)+(x0*x1)+(x1*x1)));
                    preI = i;
                end
            end
            iy = iy / 12;
            ix = ix / 12;
        end
        
    end
end

