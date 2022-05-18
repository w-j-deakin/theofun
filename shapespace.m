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
        
        %% functional landscapes

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
        
        % smoa landscape
        function [ixLS, iyLS] = SMOALandscape(obj,n)
            ix = zeros(obj.Nx, obj.Ny);
            iy = zeros(obj.Nx, obj.Ny);
            for i = 1:obj.Nx
                for j = 1:obj.Ny
                    I = (i-1)*obj.Ny + j;
                    [ix(i,j), iy(i,j)] = obj.shapes(I).SMOA(n);
                end
            end

            [xgrid,ygrid] = meshgrid(obj.X, obj.Y);
            ixLS = landscape(xgrid', ygrid', ix);
            iyLS = landscape(xgrid', ygrid', iy);
        end
    end
end

