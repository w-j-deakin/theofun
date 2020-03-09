classdef landscape
    
    % landscape class: each landscape object contains x,y and z data of
    % points in a grid, and has multiple functions for visualisation of
    % this data.
    
    properties
        nx
        ny
        nv
        xgrid
        ygrid
        zgrid
        xvec
        yvec
        zvec
        v2g
        cbarR
    end
    
    methods
        %% CONCTRUCTOR
        function obj = landscape(xgrid,ygrid,zgrid,cR)
            if nargin < 4
                cR = [min(min(zgrid)) max(max(zgrid))];
            end
            obj.cbarR = cR;
            obj.xgrid = xgrid;
            obj.ygrid = ygrid;
            obj.zgrid = zgrid;
            [obj.nx,obj.ny] = size(xgrid);
            obj.v2g = 1:obj.nx*obj.ny;
            obj.xvec = obj.vector(xgrid);
            obj.yvec = obj.vector(ygrid);
            obj.zvec = obj.vector(zgrid);
            
            k = 0;
            
            for i = 1:length(obj.xvec)
                if isnan(obj.zvec(i - k))
                    obj.xvec(i - k) = [];
                    obj.yvec(i - k) = [];
                    obj.zvec(i - k) = [];
                    obj.v2g(i - k) = [];
                    
                    k = k + 1;
                end
            end
            
            obj.nv = length(obj.xvec);
        end
        
        %% INDEXING
        
        % Get grid position from vector position
        
        function [i,j] = vec2grid(obj,I,isTrunc)
            
            if nargin < 3
                isTrunc = true;
            end
            
            nY = obj.ny;
            
            if isTrunc
                I = obj.v2g(I);
            end
            
            j = mod(I,nY);
            if j == 0
                j = nY;
            end
            k = I - j;
            i = (k/nY) + 1;
        end
        
        % Get vector position from grid position
        
        function i = grid2vec(obj,I,J,isTrunc)
            
            if nargin < 4
                isTrunc = true;
            end
            
            nY = obj.ny;
            
            k = ((I-1)*nY + J);
            
            if isTrunc
                search = true;
                i = 0;
                while search
                    i = i+1;
                    if obj.v2g(i) == k
                        search = false;
                    end
                end 
            else
                i = k;
            end
            
        end
        
        % Convert input grid to vector
        
        function v = vector(obj,grid,nX,nY)
            if nargin < 3
                nX = obj.nx;
            end
            if nargin < 4
                nY = obj.ny;
            end
            
            v = zeros(nX*nY,1);
            
            for i = 1:nX
                for j = 1:nY
                    I = obj.grid2vec(i,j);
                    v(I) = grid(i,j);
                end
            end
        end
        
        % Convert input vector to grid
        
        function g = grid(obj,v)
            
            n = length(v);
            g = nan(obj.nx,obj.ny);
            
            for i = 1:n
                [I,J] = obj.vec2grid(i);
                g(I,J) = v(i);
            end
            
        end
        
        % Delete indexed values
        
        function LS = delIDX(obj,D)
            xG = obj.xgrid;
            yG = obj.ygrid;
            zG = obj.zgrid;
            
            for i = 1:length(D)
                [I,J] = obj.vec2grid(D(i));
                zG(I,J) = nan;
            end
            
            LS = landscape(xG,yG,zG);
        end
        
        %% INTERPOLATION
        
        % Get an interpolation function using 'scatteredInterpolant'
        function z = interpolate(obj)
            z = scatteredInterpolant(obj.xvec,obj.yvec,obj.zvec);
        end
        
        % Regrid the surface to fit a new data density. (used mainly for
        % visualisation purposes).
        % 0 < factor, if factor < 1, decrease density, if > 1, increase
        % density
        function [X,Y,Z] = regrid(obj,factor)
            mx = min(min(obj.xgrid));
            Mx = max(max(obj.xgrid));
            my = min(min(obj.ygrid));
            My = max(max(obj.ygrid));
               
            nX = round(obj.nx + (obj.nx - 1)*factor);
            nY = round(obj.ny + (obj.ny - 1)*factor);
            
            [X,Y] = meshgrid(linspace(mx,Mx,nX),linspace(my,My,nY));
            Z1 = interp2(obj.xgrid',obj.ygrid',obj.zgrid',X,Y,'linear');
            Z = interp2(obj.xgrid',obj.ygrid',obj.zgrid',X,Y,'cubic');
            for i = 1:nX
                for j = 1:nY
                    if isnan(Z(j,i))
                        Z(j,i) = Z1(j,i);
                    end
                end
            end
        end
        
        % Use cubic and linear interpolation to find the surface height at
        % any x,y point.
        function z = getZ(obj,x,y)
            n = length(x);
            z = zeros(n,1);
            
            for i = 1:n
                z(i) = interp2(obj.xgrid',obj.ygrid',obj.zgrid',x(i),y(i),'cubic');
                if isnan(z(i))
                    z(i) = interp2(obj.xgrid',obj.ygrid',obj.zgrid',x(i),y(i),'linear');
                end
            end
        end
        
        % Goodness of fit of second order polynomial surface
        function gof = poly22fit(obj)
            [fitObject,gof] = fit(horzcat(obj.xvec,obj.yvec),obj.zvec,'poly22');
            figure
            plot(fitObject);
            hold on
            plot3(obj.xvec,obj.yvec,obj.zvec,'ko','MarkerFaceColor','r');
            disp(gof);
        end
        
        %% SAMPLING
        
        % Increase density
        
        function HDLS = HD(obj,nx,ny)
            
            minx = min(min(obj.xgrid));
            miny = min(min(obj.ygrid));
            maxx = max(max(obj.xgrid));
            maxy = max(max(obj.ygrid));
            
            X = linspace(minx,maxx,nx);
            Y = linspace(miny,maxy,ny);
            
            xG = zeros(nx,ny);
            yG = zeros(nx,ny);
            zG = zeros(nx,ny);
            
            Z = obj.interpolate;
            
            H = ones(obj.nx,obj.ny);
            xyz = cat(3,obj.xgrid,obj.ygrid,H);
            ptCld = pointCloud(xyz);
            
            for i = 1:nx
                for j = 1:ny
                    xG(i,j) = X(i);
                    yG(i,j) = Y(j);
                    zG(i,j) = Z(X(i),Y(j));
                    for i1 = 2:obj.nx
                        for j1 = 2:obj.ny
                            if X(i) > obj.xgrid(i1-1,j1-1) && Y(j) > obj.ygrid(i1-1,j1-1)
                                if X(i) == obj.xgrid(i1,j1) && Y(j) == obj.ygrid(i1,j1)
                                    if isnan(obj.zgrid(i1,j1))
                                        zG(i,j) = nan;
                                    end
                                elseif X(i) == obj.xgrid(i1,j1) && Y(j) < obj.ygrid(i1,j1)
                                    if isnan(obj.zgrid(i1,j1)) || isnan(obj.zgrid(i1,j1-1))
                                        zG(i,j) = nan;
                                    end
                                elseif X(i) < obj.xgrid(i1,j1) && Y(j) == obj.ygrid(i1,j1)
                                    if isnan(obj.zgrid(i1,j1)) || isnan(obj.zgrid(i1-1,j1))
                                        zG(i,j) = nan;
                                    end
                                end
                            end
                        end
                    end
                    if ~isnan(zG(i,j))
                        I = findNearestNeighbors(ptCld,[X(i) Y(j) 1], 4);
                        for k = 1:4
                            if isnan(obj.zgrid(I(k)))
                                zG(i,j) = nan;
                            end
                        end
                    end
                end
            end
            
            HDLS = landscape(xG,yG,zG);
        end
        
        % Bootstrap surface
        
        function [LS,ri,r] = bootstrap(obj,ri)
            
            n = obj.nv;
            nX = obj.nx;
            nY = obj.ny;
            if nargin < 2
                ri = randi(n,nX,nY);
            end
            xG = zeros(nX,nY);
            yG = zeros(nX,nY);
            zG = zeros(nX,nY);
            
            for i = 1:nX
                for j = 1:nY
                    if isnan(obj.zgrid(i,j))
                        zG(i,j) = nan;
                        xG(i,j) = nan;
                        yG(i,j) = nan;
                    else
                        zG(i,j) = obj.zvec(ri(i,j));
                        xG(i,j) = obj.xvec(ri(i,j));
                        xG(i,j) = obj.yvec(ri(i,j));
                    end
                end
            end
            LS = landscape(xG,yG,zG);
            r = LS.vector(ri);
        end

        
        %% PLOTS
        
        % Get x,y difference between regular intervals on the grid
        function [rx,ry] = diff(obj)
            rx = abs(obj.xgrid(2,1) - obj.xgrid(1,1));
            ry = abs(obj.ygrid(1,2) - obj.ygrid(1,1));
        end
        
        % Plot 1 individual pixel colour corresponding to z value 
        
        %(NOT RECOMMENDED)
        function colourSqr(obj,x,y)
            if isnan(obj.zgrid(x,y))
                return
            else
                ctr = [obj.xgrid(x,y) obj.ygrid(x,y)];
                [rx,ry] = obj.diff;
                LX = ctr(1) - rx/2;
                RX = ctr(1) + rx/2;
                DY = ctr(2) - ry/2;
                UY = ctr(2) + ry/2;
                c = colormap;
                rc = obj.cbarR(2) - obj.cbarR(1);
                cP = (obj.zgrid(x,y) - obj.cbarR(1))/rc;
                [nC,~] = size(c);
                cI = round(cP*(nC-1))+1;
                col = c(cI,:);
                fill([LX LX RX RX],[DY UY UY DY],col,'LineStyle','none');
            end
        end
        
        % Plot Heatmap as a rasterised image in the axes space
        function fillHM(obj)
            mx = min(min(obj.xgrid));
            Mx = max(max(obj.xgrid));
            my = min(min(obj.ygrid));
            My = max(max(obj.ygrid));
            lim = [mx Mx;
                   my My];
           
            [X,Y,Z] = obj.regrid(1);
            
            [nY,nX] = size(Z);
            
            C = zeros(nY,nX,3);
            cmap = colormap;
            [nc,~] = size(cmap);
            rc = obj.cbarR(2) - obj.cbarR(1);
            z = Z-obj.cbarR(1);
            z = z./rc;
            
            colZ = linspace(0,1,nc);
            
            for i = 1:nX
                for j = 1:nY
                    if isnan(z(j,i))
                        C(j,i,:) = [0.75 0.75 0.75];
                    elseif z(j,i) >= 1
                        C(j,i,:) = cmap(end,:);
                    elseif z(j,i) <= 0
                        C(j,i,:) = cmap(1,:);
                    else
                        low = ceil(z(j,i)*(nc-1));
                        high = low+1;
                        lowC = cmap(low,:);
                        highC = cmap(high,:);
                        r = highC - lowC;
                        rz = colZ(high) - colZ(low);
                        cz = (z(j,i) - colZ(low))/rz;
                        R = r*cz;

                        C(j,i,:) = lowC + R;
                    end
                end
            end
            
            imagesc(lim(1,:),lim(2,:),C);
            set(gca,'YDir','normal');
        end
        
        % Plot Contourmap
        function fillCM(obj,cont,zlim)
            if nargin < 3
                zlim = [-inf inf];
            end
            
            mx = min(min(obj.xgrid));
            Mx = max(max(obj.xgrid));
            my = min(min(obj.ygrid));
            My = max(max(obj.ygrid));

            [X,Y,Z] = obj.regrid(100);
            
            [nY,nX] = size(Z);
            z = ones(nY,nX);
            z = z*obj.cbarR(2);
            
            for i = 1:nY
                for j = 1:nX
                    if Z(i,j) < zlim(1)
                        Z(i,j) = zlim(1);
                    elseif Z(i,j) > zlim(2)
                        Z(i,j) = zlim(2);
                    end
                    if isnan(Z(i,j))
                        z(i,j) = obj.cbarR(1);
                    end
                end
            end
            
            contourf(X,Y,z,[obj.cbarR(1) obj.cbarR(2)],'LineWidth',1.1);
            hold on
            fill([mx mx Mx Mx],[my My My my],'k','LineStyle','none','FaceAlpha',0.2);
            contourf(X,Y,Z,cont);
        end
        
        % Function to create a contourmap figure
        
        function heatMap(obj,hue,numCont,ax,pcx,pcy,isSmall)
            
            if nargin < 7
                isSmall = false;
            end
            if nargin < 6
                pcy = 2;
            end
            if nargin < 5
                pcx = 1;
            end
            
            mx = min(min(obj.xgrid));
            Mx = max(max(obj.xgrid));
            my = min(min(obj.ygrid));
            My = max(max(obj.ygrid));
            
            xL = [mx Mx];
            yL = [my My];
            
            X = [xL(1) xL(1) xL(2) xL(2)];
            Y = [yL(1) yL(2) yL(2) yL(1)];
            
            if nargin < 4
                ax = axes;
            end
            position = ax.Position;
            
            mx = ceil(mx*100)/100;
            Mx = floor(Mx*100)/100;
            my = ceil(my*100)/100;
            My = floor(My*100)/100;
            
            cX = position(3)/20;
            
                
            ax.Position = [position(1) position(2) position(3)-2*cX position(4)];
            axis equal
            cP = [position(1)+19*cX position(2) cX position(4)];
            cont = linspace(min(obj.zvec),max(obj.zvec),numCont+1);
            map = obj.HueColorMap(hue);
            colormap(map);
            caxis(obj.cbarR);
            obj.fillCM(cont(1:numCont),[min(obj.zvec) max(obj.zvec)]);
            
            hold on
            fill(X,Y,'w','FaceAlpha',0);
            xlim(xL);
            ylim(yL);

%             cAx = axes;
%             cAx.Position = cP;
%             range = [min(min(obj.zgrid)) max(max(obj.zgrid))];
%             obj.DrawCBar(cAx,range,nan,isSmall);
%             axes(ax);
%             
            if isSmall
                axis off
            else
                set(ax,'xtick',[mx (mx+Mx)/2 Mx]);
                S1 = ['PC',num2str(pcx)];
                set(ax,'xticklabel', {num2str(mx), S1, num2str(Mx)});
                set(ax,'ytick',[my (my+My)/2 My]);
                S2 = ['PC',num2str(pcy)];
                set(ax,'yticklabel',{num2str(my), S2, num2str(My)});
            end
        end
        
        
        %% COLORMAPS
        
        % Colorblind friendly colormap, single color specified by input
        % hue.
        function map = HueColorMap(obj,hue)
            map = zeros(256,3);
            
            sat = linspace(0,1,256);
            val = linspace(1,0.75,256);
            
            for i = 1:256
                HSV = [hue sat(i) val(i)];
                RGB = hsv2rgb(HSV);
                map(i,:) = RGB;
            end
        end
        
        % Test a colormap
        function testColorMap(obj,hue)
            map = obj.HueColorMap(hue);
            figure
            X = [0 0 1 1];
            for i = 1:10
                fill(X,[i+0.1 i+0.9 i+0.9 i+0.1],map(i*200,:));
                hold on
            end
        end
        
        % Draw a colorbar (with max-min)
        function DrawCBar(obj,ax,range,bars,isSmall)
            if nargin < 5
                isSmall = false;
            end
            axes(ax);
            cmap = colormap;
            [N,~] = size(cmap);
            
            Y = linspace(obj.cbarR(1),obj.cbarR(2),N);
            d = Y(2) - Y(1);
            Y = Y - d/2;
            Y = horzcat(Y,Y(end)+d/2);
            x = [0 0 1 1];
            fill(x,[Y(1) Y(end) Y(end) Y(1)],'w','LineStyle','none');
            hold on
            if nargin < 3
                range = [Y(1) Y(end)];
            end
            if nargin < 4
                bars = nan;
            end
            k = 0;
            for i = 1:N
                y = [Y(i) Y(i+1) Y(i+1) Y(i)];
                if Y(i+1) < range(1) || Y(i) > range(2)
                    fill(x,y,cmap(i,:),'LineStyle','none','FaceAlpha',0.25);
                else
                    fill(x,y,cmap(i,:),'LineStyle','none');
                    k = k+1;
                    yL(k) = i;
                end
                hold on
            end
            fill(x,[Y(yL(1)) Y(yL(end)+1) Y(yL(end)+1) Y(yL(1))],'k','FaceAlpha',0);
            for i = 1:length(bars)
                if isSmall
                    plot([0 1],[bars(i) bars(i)],'k-','LineWidth',1.1);
                    hold on
                    plot([0.1 0.9],[bars(i) bars(i)],'w-','LineWidth',0.6);
                else
                    plot([0 1],[bars(i) bars(i)],'k-','LineWidth',2);
                    hold on
                    plot([0.1 0.9],[bars(i) bars(i)],'w-','LineWidth',1);
                end
            end
            xlim([0 1]);
            ylim([Y(1) Y(end)]);
            ax.YAxisLocation = 'right';
            ax.XColor = 'none';
            ax.Color = 'none';
            ax.Box = 'off';
            
        end
        
        function v = GetVector(obj,x,y)
            for i = 1:obj.nx-1
                if x >= obj.xgrid(i,1) && x <= obj.xgrid(i+1,1)
                    X = i;
                end
            end
            for i = 1:obj.ny-1
                if y >= obj.ygrid(1,i) && y <= obj.ygrid(1,i+1)
                    Y = i;
                end
            end
            D = [-1 0 1];
            V = zeros(9,2);
            if isnan(obj.zgrid(X,Y))
                v = nan;
            else
                for i = 1:3
                    for j = 1:3
                        idx = (i-1)*3 + j;
                        if isnan(obj.zgrid(X+D(i),Y+D(j)))
                            v = nan;
                            return
                        else
                            M = obj.zgrid(X+D(i),Y+D(j)) - obj.zgrid(X,Y);
                        end
                        V(idx,:) = [D(i) D(j)] * M;
                    end
                end
            end
            v = sum(V);
        end

        %%
    end
end

