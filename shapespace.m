classdef shapespace
    
    % Shape Space built from a PC morphospace. Contains a grid of
    % theoretical morphologies.
    
    properties
        taxa
        % taxaSet shapespace is built from
        shapes = theoShapeN.empty;
        % theoretical morphologies
        pcx
        % x axis PC
        pcy
        % y axis PC
        Nx
        % number of shapes drawn along x axis
        Ny
        % number of shapes drawn along x axis
        xR
        % x range ([minx maxx])
        yR
        % y range ([miny maxy])
        X
        % X vector of grid positions
        Y
        % Y vector of grid positions
    end
    
    methods
        %% CONSTRUCTOR
        function obj = shapespace(taxaSet,PCx,PCy,n,xR,yR)
            obj.taxa = taxaSet;
            obj.pcx = PCx;
            obj.pcy = PCy;
            obj.xR = xR;
            m = taxaSet.mu;
            
            rx = xR(2) - xR(1);
            ry = yR(2) - yR(1);
            
            r = rx/ry;
            
            [nr,nc] = size(n);
            if nr == 1 && nc == 1
                ny = round(sqrt(n/r));
                nx = floor(ny*r);
            else
                nx = n(1,1);
                ny = round(nx/r);
            end
            
            R = rx/(nx-1);
            
            obj.Nx = nx;
            
            xC = zeros(nx,1);
            
            for i = 1:nx
                xC(i) = xR(1) + (i-1)*R;
            end
            
            obj.X = xC;
            
            ny = 0;
            run = true;
            while run
                ny = ny + 1;
                yC(ny) = yR(1) + (ny-1)*R;
                if yC(ny) >= yR(2)
                    run = false;
                end
            end
            over = yC(ny) - yR(2);
            over = over/2;
            yC = yC - over;
            
            obj.Y = yC;
            obj.Ny = ny;
            obj.yR = [min(yC) max(yC)];
            
            coef = vertcat(taxaSet.coeff(:,PCx)',taxaSet.coeff(:,PCy)');
            
            for i = 1:nx
                for j = 1:ny
                    H = ([xC(i) yC(j)]*coef) + m;
                    obj.shapes((i-1)*ny + j) = theoShapeN([xC(i) yC(j)],H);
                end
            end
        end
        
        %% PLOTS
        
        % plot individual background shapes
        
        function s = plotGrid(obj,col)
            if nargin < 2
                col = [0.75 0.75 0.75];
            end
            
            n = round(500/sqrt(obj.Nx*obj.Ny));
            s = (obj.X(2) - obj.X(1))/3;
            
            for i = 1:obj.Nx
                for j = 1:obj.Ny
                    I = (i-1)*obj.Ny + j;
                    [x,y] = obj.shapes(I).draw(n,obj.X(i),obj.Y(j),s);
                    fill(x,y,col,'LineStyle','none');
                    hold on
                end
            end
        end
        
        % plot ShapeSpace
        
        function Plot(obj)
            nx = obj.Nx;
            ny = obj.Ny;
            XR = obj.xR;
            YR = obj.yR;
            ax1PC = obj.pcx;
            ax2PC = obj.pcy;
            
            bX = (XR(2) - XR(1)) / (nx-1);
            bY = (YR(2) - YR(1)) / (ny-1);
            
            xL = [(XR(1) - bX) (XR(2) + bX)];
            yL = [(YR(1) - bY) (YR(2) + bY)];
            
            figure
            
            obj.taxa.Cmorph(ax1PC,ax2PC,false);
            xlim(xL);
            ylim(yL);
            
            
            ax = gca;
            p = ax.Position;
            
            pX = zeros(nx,1);
            pY = zeros(ny,1);
            
            rx = p(3)/(nx+1);
            ry = p(4)/(ny+1);
            
            pX(1) = p(1) + 0.55 * rx;
            pY(1) = p(2) + 0.55 * ry;
            
            for i = 2:nx
                pX(i) = pX(1) + (i-1) * rx;
            end
            
            for i = 2:ny
                pY(i) = pY(1) + (i-1) * ry;
            end
            
            w = rx * 0.9;
            h = ry * 0.9;
            
            obj.plotGrid;
            
            axes;
            obj.taxa.Cmorph(ax1PC,ax2PC,true);
            xlim(xL);
            ylim(yL);
            
            axis off
        end
        
        % Draw specific outline
        
        function [ox,oy] = drawOL(obj,x,y)
            
            
            idx = (x-1)*(obj.Ny) + y;
            
            [ox,oy] = obj.shapes(idx).draw(500);
            fill(ox,oy,col,'LineStyle','none');
        end
        
         %% density landscape

        function ls = DensityLandscape(obj,allAxes)
            
            data = obj.taxa.scores;
            % bandwidth using silverman's rule of thumb
            if nargin < 2
                allAxes = true;
            end

            if allAxes
                pcX = obj.pcx;
                pcY = obj.pcy;
            else
                data = horzcat(data(:,obj.pcx),data(:,obj.pcy));
                pcX = 1;
                pcY = 2;
            end

            [nT,nV] = size(data);

            bw = zeros(nV,1);
            for i = 1:nV
                sig = std(data(:,i));
                bw(i) = sig * ( 4 / ((nV + 2) * nT) ) ^ ( 1 / (nV + 4) );
            end

            nS = length(obj.shapes);
            P = zeros(nS,nV);
            for i = 1:length(obj.shapes)
                P(i,pcX) = obj.shapes(i).v(1);
                P(i,pcY) = obj.shapes(i).v(2);
            end
            d = mvksdensity(data,P,'Bandwidth',bw);
            D = zeros(obj.Nx,obj.Ny);
            for i = 1:obj.Nx
                for j = 1:obj.Ny
                    I = (i-1)*obj.Ny + j;
                    if obj.shapes(I).isIntersect
                        D(i,j) = nan;
                    else
                        D(i,j) = d(I);
                    end
                end
            end

            [xgrid,ygrid] = meshgrid(obj.X, obj.Y);
            ls = landscape(xgrid', ygrid', D);
        end
        
        %% functional landscapes

        % area landscape
        function ls = AreaLandscape(obj,n)
            if nargin < 2
                n = 500;
            end
            z = zeros(obj.Nx, obj.Ny);
            for i = 1:obj.Nx
                for j = 1:obj.Ny
                    I = (i-1)*obj.Ny + j;
                    [x,y] = obj.shapes(I).Draw(n,0,0,1);
                    z(i,j) = polyarea(x,y);
                end
            end
            
            [xgrid,ygrid] = meshgrid(obj.X, obj.Y);
            ls = landscape(xgrid', ygrid', z);
        end

        % wing aspect ratio landscape
        function ls = AspectRatioLandscape(obj,n)
            z = zeros(obj.Nx, obj.Ny);
            for i = 1:obj.Nx
                for j = 1:obj.Ny
                    I = (i-1)*obj.Ny + j;
                    z(i,j) = obj.shapes(I).AspectRatio(n);
                end
            end

            [xgrid,ygrid] = meshgrid(obj.X, obj.Y);
            ls = landscape(xgrid', ygrid', z);
        end
        
        % moa landscapes
        function [i1xLS, i1yLS, i2xLS, i2yLS, i2zLS] = MoALandscape(obj,n)
            i1x = zeros(obj.Nx, obj.Ny);
            i1y = zeros(obj.Nx, obj.Ny);
            i2x = zeros(obj.Nx, obj.Ny);
            i2y = zeros(obj.Nx, obj.Ny);
            i2z = zeros(obj.Nx, obj.Ny);

            for i = 1:obj.Nx
                for j = 1:obj.Ny
                    I = (i-1)*obj.Ny + j;
                    [i1x(i,j), i1y(i,j), i2x(i,j), i2y(i,j), i2z(i,j)] = obj.shapes(I).MoA(n);
                end
            end

            [xgrid,ygrid] = meshgrid(obj.X, obj.Y);
            i1xLS = landscape(xgrid', ygrid', i1x);
            i1yLS = landscape(xgrid', ygrid', i1y);
            i2xLS = landscape(xgrid', ygrid', i2x);
            i2yLS = landscape(xgrid', ygrid', i2y);
            i2zLS = landscape(xgrid', ygrid', i2z);
        end

        % muscle attachment landscape
        function [maMLS, ma5LS, ma95LS] = MALandscape(obj,L,alpha,nSamp,LRange,alphaRange)
            xgrid = zeros(obj.Nx, obj.Ny);
            ygrid = zeros(obj.Nx, obj.Ny);
            zgrid = zeros(obj.Nx, obj.Ny);

            for i = 1:obj.Nx
                for j = 1:obj.Ny
                    xgrid(i,j) = obj.X(i);
                    ygrid(i,j) = obj.Y(j);
                end
            end

            if nargin < 4
                for i = 1:obj.Nx
                    for j = 1:obj.Ny
                        I = (i-1)*obj.Ny + j;
                        zgrid(i,j) = obj.shapes(I).MuscleAttachmentArea(L,alpha);
                    end
                end
                maMLS = landscape(xgrid,ygrid,zgrid);
            else
                zgrid5 = zeros(obj.Nx, obj.Ny);
                zgrid95 = zeros(obj.Nx, obj.Ny);
                a = zeros(nSamp,1);
                for i = 1:obj.Nx
                    for j = 1:obj.Ny
                        I = (i-1)*obj.Ny + j;
                        for k = 1:nSamp
                            alp = alpha + (rand-0.5) * alphaRange;
                            len = L + (rand-0.5) * LRange;
                            a(k) = obj.shapes(I).MuscleAttachmentArea(len,alp);
                        end
                        zgrid(i,j) = nanmean(a);
                        zgrid5(i,j) = prctile(a,5);
                        zgrid95(i,j) = prctile(a,95);
                    end
                end
                maMLS = landscape(xgrid,ygrid,zgrid);
                ma5LS = landscape(xgrid,ygrid,zgrid5);
                ma95LS = landscape(xgrid,ygrid,zgrid95);
            end
        end

        % wing agility landscape (Harvey et al 2022)
        function dqLS = AgilityLandscape(obj,n)
            dq = zeros(obj.Nx,obj.Ny);

            for i = 1:obj.Nx
                for j = 1:obj.Ny
                    I = (i-1)*obj.Ny + j;
                    dq(i,j) = obj.shapes(I).WingAgility(n);
                end
            end

            [xgrid,ygrid] = meshgrid(obj.X, obj.Y);
            dqLS = landscape(xgrid', ygrid', dq);
            
        end
    end
end

