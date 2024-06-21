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

        function [i1x, i1y, i2x, i2y, i2z, i3y] = MoA(obj, n)
            if (obj.isIntersect)
                i1x = nan;
                i1y = nan;
                i2x = nan;
                i2y = nan;
                i2z = nan;
                i3y = nan;
            else
                [x,y] = obj.draw(n,0,0,1);
                [~,baseID] = min(x);
                [~,tipID] = max(x);
                [~,a] = findAngle2D([1 0], [(x(tipID) - x(baseID)) (y(tipID) - y(baseID))]);

                if y(tipID)-y(baseID) > 0
                    a = -a;
                end
                x = x - x(baseID);
                y = y - y(baseID);

                points = Rotate([x y], a);

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
                    i3y = i3y + A * ((x0 + x1)^3 + (x0 + x1)*(x1-x0)^2);
                    preI = i;
                end
                
                i1x = abs(i1x / 6);
                i1y = abs(i1y / 6);
                i2x = abs(i2x / 12);
                i2y = abs(i2y / 12);
                i3y = abs(i3y / 40);
                i2z = i2x + i2y;
            end
        end

        % Scaled MoA
        function [i1y,i2y,i3y] = ScaledMoA(obj, n)
            [~,i1y,~i2y,~,i3y] = obj.MoA(n);
            [x, y] = obj.draw(n, 0, 0, 1);
            [~, baseID] = min(x);
            [~, tipID] = max(x);
            [~, a] = findAngle2D([1 0], [(x(tipID) - x(baseID)) (y(tipID) - y(baseID))]);

            if y(tipID) - y(baseID) > 0
                a = -a;
            end
            x = x - x(baseID);
            y = y - y(baseID);
            points = Rotate([x y], a);
            A = polyarea(points(:, 1), points(:, 2));
            points = points / sqrt(A);
            L = norm(max(points(:, 1)) - min(points(:, 1)));
            i1y = (i1y / L);
            i2y = (i2y / (L^2))^0.5;
            i3y = (i3y / (L^3))^(1/3);
        end

        % Muscle attachment
        function A = MuscleAttachmentArea(obj,L,alpha)
            % Cut the shape in half with a line that has an angle alpha
            % from the vertical, and intersects the top surface of the
            % outline a distance L away from the minimum x point.

            if obj.isIntersect
                A = nan;
                return
            end

            [x,y] = obj.draw(500,0,0,1);
            shp = polyshape(x,y);
            totA = area(shp);

            x = x ./ sqrt(totA);
            y = y ./ sqrt(totA);
            
            minx = min(x);
            maxx = max(x);
            r = maxx - minx;
            L = r * L;
            xP = minx + L;
            points = [0 0];

            for i = 1:499
                if (x(i) >= xP && x(i+1) < xP) || (x(i) < xP && x(i+1) >= xP)
                    p = [x(i+1)-x(i) y(i+1)-y(i)];
                    p = p / (x(i+1) - x(i));
                    p = p * (xP-x(i));
                    p = [x(i) y(i)] + p;
                    points = vertcat(points,p);
                end
            end
            points(1,:) = [];
            if isempty(points)
                A = 1;
            else
                [~,id] = max(points(:,2));
                P = points(id,:);
                x = x - P(1);
                y = y - P(2);
                ps = Rotate([x y],alpha*-1);
                ps1 = ps;
                for i = 500:-1:1
                    if ps1(i,1) > 0
                        ps1(i,:) = [];
                    end
                end
                shp1 = polyshape(ps1(:,1),ps1(:,2));
                A = area(shp1);
           
            end
        end
        
        % Wing Agility
        function dq = WingAgility(obj,res)
            if obj.isIntersect
                dq = nan;
                return;
            end
            % draw outline
            [x,y] = obj.draw(res,0,0,1);

            % find minimum x point, set to [0 0].
            [~,baseID] = min(x);
            [~,tipID] = max(x);
            [~,a] = findAngle2D([1 0], [(x(tipID) - x(baseID)) (y(tipID) - y(baseID))]);

            if y(tipID)-y(baseID) > 0
                a = -a;
            end

            x = x - x(baseID);
            y = y - y(baseID);

            points = Rotate([x y], a);

            A = polyarea(points(:,1),points(:,2));
            points = points / sqrt(A);

            x = points(:,1);
            y = points(:,2);

            % evenly space sample points along wing length
            X = linspace(0,max(x),res+2);
            dx = X(2)-X(1);
            X(1) = [];
            X(end) = [];
            Y = [0 0];

            c = zeros(res,1);

            nomYC4 = 0;
            denomYC4 = 0;

            % for each sample point
            for i = 1:res
                % for each outline segment
                ylist = 0;
                for j = 1:res
                    % find end point of segment
                    next = mod(j,res) + 1;

                    % if the current sample point lies between the x values
                    % of the segment
                    if (x(j) >= X(i) && x(next) <= X(i)) || (x(j) <= X(i) && x(next) >= X(i))
                        L = [x(j) y(j)] - [x(next) y(next)];
                        d = X(i) - x(next);
                        r = d/L(1);
                        L = L * r;
                        L = L + [x(next) y(next)];
                        ylist = vertcat(ylist,L(2));
                    end
                end
                ylist(1) = [];
                Y(1) = min(ylist);
                Y(2) = max(ylist);



                c(i) = Y(2) - Y(1);
                yc4 = Y(2) - (c(i)/4);

                yc4 = abs(yc4);
                
                nomYC4 = nomYC4 + c(i) * yc4 * dx;
                denomYC4 = denomYC4 + c(i) * dx;
            end

            yc4 = nomYC4 / denomYC4;
            crmax = max(c);

            shape = polyshape(x,y);
            [~,CGy] = centroid(shape);

            CGy = abs(CGy);

            [~,~,I2x] = MoA(obj,res);



            % formula (Harvey et al., 2022)
            dq = ((((yc4/crmax)^0.8)*crmax)-CGy) / I2x;
        end
    end
end

