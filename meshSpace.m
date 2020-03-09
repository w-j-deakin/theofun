classdef meshSpace
    
    % Meshspace, based on shapespace SS. Each shape in SS is automatically
    % meshed to an approximate number of elements, which can then be
    % functionally tested.
    
    properties
        SS
        % input shapespace
        mesh
        % array of 2D mesh objects
        xgrid
        % grid of x co-ordinates
        ygrid
        % grid of y co-ordinates
        isMeshed
        % array of boolean variables - shows which shapes in the shapespace
        % could or could not be meshed (0 = not meshed, 1 = meshed);
    end
    
    methods
        %% CONSTRUCTOR
        function obj = meshSpace(SS,nEle)
            obj.SS = SS;
            
            xgrid = zeros(SS.Nx,SS.Ny);
            ygrid = zeros(SS.Nx,SS.Ny);
            obj.isMeshed = zeros(SS.Nx, SS.Ny);
            
            for i = 1:SS.Nx
                for j = 1:SS.Ny
                    xgrid(i,j) = SS.X(i);
                    ygrid(i,j) = SS.Y(j);
                    
                    k = (i-1)*SS.Ny + j;
                    
                    if ~(SS.shapes(k).isIntersect)
                        obj.isMeshed(i,j) = 1;
                    end
                end
            end
            
            obj.mesh = meshList(SS.shapes,nEle);
            
            obj.xgrid = xgrid;
            obj.ygrid = ygrid;
        end
        
        %% INDEXING
        
        % get grid position form vector
        
        function [i,j] = vec2grid(obj,I)
            ny = obj.SS.Ny;
            
            j = mod(I,ny);
            if j == 0
                j = ny;
            end
            k = I - j;
            i = (k/ny) + 1;
        end
        
        % get vector position from grid
        
        function i = grid2vec(obj,I,J)
            ny = obj.SS.Ny;
            
            i = ((I-1)*ny + J);
        end
        
        % convert vector to grid
        
        function g = grid(obj,v)
            g = zeros(obj.SS.Nx,obj.SS.Ny);
            
            for i = 1:obj.SS.Nx
                for j = 1:obj.SS.Ny
                    I = obj.grid2vec(i,j);
                    g(i,j) = v(I);
                end
            end
        end
        
        %% FUNCTIONAL CALCULATIONS
        
        % FEA
        
        function [vms] = autoFEA(obj,t,E,v)
            tic
            if nargin < 2
                t = 0.01;
            end
            if nargin < 3
                E = 2*(10^9);
            end
            if nargin < 4
                v = 0.3;
            end
            
            nM = obj.SS.Nx*obj.SS.Ny;
            vms = obj.mesh.autoFEA(t,E,v);
            
            for i = 1:nM
                [I,J] = obj.vec2grid(i);
                if obj.isMeshed(I,J)
                    if i == 1
                        vms = vertcat(nan,vms);
                    else
                        vms = vertcat(vms(1:i-1),nan,vms(i:end));
                    end 
                end
            end
            vms = vms';
            toc
        end
        
        % FEA with random conditions
        
        function [vms] = randFEA(obj,n,r,aR,t,E,v)
            tic
            if nargin < 3
                r = 0.1;
            end
            if nargin < 4
                aR = [-pi/4 pi/4];
            end
            if nargin < 5
                t = 0.01;
            end
            if nargin < 6
                E = 2*(10^9);
            end
            if nargin < 7
                v = 0.3;
            end
            
            nM = obj.SS.Nx*obj.SS.Ny;
            vms = obj.mesh.randFEA(n,r,aR,t,E,v);
            
            for i = 1:nM
                [I,J] = obj.vec2grid(i);
                if ~obj.isMeshed(I,J)
                    if i == 1
                        vms = horzcat(nan(n,1),vms);
                    else
                        vms = horzcat(vms(:,1:i-1),nan(n,1),vms(:,i:end));
                    end 
                end
            end
            
%             writematrix(vms,'vmsBtstrp.csv');
            
            toc
        end
        
        % rotational inertia and tip velocity
        
        function [v,RI] = autoRI(obj)
            tic
            
            nM = obj.SS.Nx*obj.SS.Ny;
            [RI,v] = obj.mesh.autoRI;
            
            for i = 1:nM
                [I,J] = obj.vec2grid(i);
                if obj.isMeshed(I,J)
                    if i == 1
                        v = vertcat(nan,v);
                        RI = vertcat(nan,RI);
                    else
                        v = vertcat(v(1:i-1),nan,v(i:end));
                        RI = vertcat(RI(1:i-1),nan,RI(i:end));
                    end 
                end
            end
            v = v';
            RI = RI';
            toc
        end
        
        % RI and v with random conditions
        
        function [v,RI] = randRI(obj,n)
            tic
            nM = obj.SS.Nx*obj.SS.Ny;
            [RI,v] = obj.mesh.randRI(n);
            
            for i = 1:nM
                [I,J] = obj.vec2grid(i);
                if ~obj.isMeshed(I,J)
                    if i == 1
                        v = horzcat(nan(n,1),v);
                        RI = horzcat(nan(n,1),RI);
                    else
                        v = horzcat(v(:,1:i-1),nan(n,1),v(:,i:end));
                        RI = horzcat(RI(:,1:i-1),nan(n,1),RI(:,i:end));
                    end 
                end
            end
            toc
        end
        
        %% PLOTS
        
        % generate landscape
        
        function LS = saveLS(obj,zgrid)
            LS = landscape(obj.xgrid,obj.ygrid,zgrid);
        end
        
        % generate median, 95 CL and 5 CL landscapes of random data
        
        function [LSm,LS5,LS95] = saveRandLS(obj,var)
            
            vm = nanmean(var);
            vp = prctile(var,[5 95],1);
            
            v5 = vp(1,:);
            v95 = vp(2,:);
            
            LSm = landscape(obj.xgrid,obj.ygrid,obj.grid(vm),[min(v5) max(v95)]);
            LS5 = landscape(obj.xgrid,obj.ygrid,obj.grid(v5),[min(v5) max(v95)]);
            LS95 = landscape(obj.xgrid,obj.ygrid,obj.grid(v95),[min(v5) max(v95)]);
            
        end
        
        % FEA landscape
        
        function LS = vmsHeat(obj)
             vms = obj.autoFEA;
             Gvms = obj.grid(vms);
            
             LS = obj.saveLS(Gvms);
        end
        
        % tip velocity landscape
        
        function LS = rotHeat(obj)
            v = obj.autoRI;
            Gv = obj.grid(v);

            LS = obj.saveLS(Gv);
        end
        
        % pareto front of two landscapes, with taxon scores superimposed
        
        function taxaPF(obj,LS1,LS2,opt1,opt2,x,y)
            
            if nargin < 6
                x = obj.SS.taxa.scores(:,obj.SS.pcx);
            end
            if nargin < 7
                y = obj.SS.taxa.scores(:,obj.SS.pcy);
            end
            if nargin < 4
                opt1 = true;
            end
            if nargin < 5
                opt2 = true;
            end
            
            z1 = LS1.getZ(x,y);
            z2 = LS2.getZ(x,y);
            
            r1 = max(z1) - min(z1);
            r2 = max(z2) - min(z2);
            
            b1 = r1*0.1;
            b2 = r2*0.1;
            
%             ax = axes; 
            
            p = pareto(LS1.zvec,LS2.zvec,opt1,opt2);
            pT = pareto(z1,z2,opt1,opt2);
            p.Plot;
            hold on
            pT.Plot([1 0 0],true);
            xL = [min(z1)-b1 max(z1)+b1];
            yL = [min(z2)-b2 max(z2)+b2];
            plot([xL(1) xL(1) xL(2) xL(2) xL(1)],[yL(1) yL(2) yL(2) yL(1) yL(1)],'r-');
            xlabel('VMS');
            ylabel('RE');
        end
        
        % Pareto optimal landscape
        
        function LS = morphPF(obj,LS1,LS2,opt1,opt2,hue,showTaxa)
            
            if nargin < 7
                showTaxa = true;
            end
            
            
            p = pareto(LS1.zvec,LS2.zvec,opt1,opt2);
            R = p.RankRatio;
            
            LS = landscape(LS1.xgrid,LS1.ygrid,LS1.grid(R));
            LS.heatMap(hue);
            axis equal
            
            if showTaxa
                hold on
                obj.SS.taxa.Plot(obj.SS.pcx,obj.SS.pcy,'ko','w',false);
            end
            
        end
        
        
        function PlotFEM(obj)
            figure
            k = 0;
            s = obj.SS.plotGrid;
            s = s*0.75;
            close all;
            figure
            for i = 1:obj.SS.Nx
                for j = 1:obj.SS.Ny
                    if obj.isMeshed(i,j)
                        k = k+1;
                        obj.mesh.mesh(k).DrawFEM([obj.xgrid(i,j) obj.ygrid(i,j)],s);
                        hold on
                    end
                    disp((i-1)*obj.SS.Ny + j);
                end
            end
            axis equal
        end
            
        %%
    end
end

