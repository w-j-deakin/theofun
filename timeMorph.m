classdef timeMorph
    
    % Time Split Morphospace object
    
    
    properties
        nBins
        bins
        allTaxa
        binTaxa = taxaSet.empty;
        binName
        binCol
    end
    
    methods
        %% CONSTRUCTOR
        function obj = timeMorph(taxa,Times,Names,Colours)
            timeSplit = taxa.timeRange;
            
            times = importdata(Times);
            
            names = importdata(Names);
            
            colours = importdata(Colours);
            
            minT = timeSplit(1);
            maxT = timeSplit(2);
            
            obj.nBins = maxT - minT + 1;
            obj.allTaxa = taxa;
            
            if nargin < 2
                obj.bins = 0:1:obj.nBins;
            else
                [rT,cT] = size(times);
                if rT > 1 && cT > 1
                    error('Input times must be a vector');
                elseif rT < (obj.nBins+1) && cT < (obj.nBins+1)
                    error('Input time vector must have more elements than number of bins');
                else
                    obj.bins = times(1:obj.nBins+1);
                end
            end
            
            if nargin < 3
                cellNames = cell(obj.nBins,1);
                
                for i = 1:obj.nBins
                    cellNames{i} = num2str(i);
                end
                
                obj.binName = string(cellNames);
            else
                [rN,cN] = size(names);
                if rN > 1 && cN > 1
                    error('Name vector must have only 1 row or column');
                elseif rN < obj.nBins && cN < obj.nBins
                    error('Name vector must be equal to or longer than number of bins');
                else
                    for i = 1:obj.nBins
                        obj.binName{i} = names{i};
                    end
                end
            end
            
            if nargin < 4
                obj.binCol = zeros(obj.nBins,3);
                
                colormap hsv
                cmap = colormap;
                [nC,~] = size(cmap);
                
                close all
                r = nC/obj.nBins;
                
                for i = 1:obj.nBins
                    obj.binCol(i,:) = cmap(i*r,:);
                end
            else
                [nC,cC] = size(colours);
                if nC < obj.nBins
                    error('Colour array must have at least as many rows as there are time bins');
                elseif cC ~= 3
                    error('Colour array must have 3 columns (RGB)');
                else
                    obj.binCol = colours(1:obj.nBins,:);
                end
            end
            
            for i = 1:obj.nBins
                obj.binTaxa(i) = taxa.timeSplit(minT + i - 1);
            end
        end
        
        %% PLOTS
        function cladeDraw(obj,xPC,yPC,isBckgrnd)
            figure
            
            if nargin < 4
               isBckgrnd = true; 
            end
            
            taxa = obj.allTaxa;
            
            minX = min(taxa.scores(:,xPC));
            minY = min(taxa.scores(:,yPC));
            maxX = max(taxa.scores(:,xPC));
            maxY = max(taxa.scores(:,yPC));
            
            xB = (maxX - minX) / 20;
            yB = (maxY - minY) / 20;
            
            xL = [(minX - xB) (maxX + xB)];
            yL = [(minY - yB) (maxY + yB)];
            
            X = [xL(1) xL(1) xL(2) xL(2)];
            Y = [yL(1) yL(2) yL(2) yL(1)];
            
            if isprime(obj.nBins)
                nP = obj.nBins + 1;
            else
                nP = obj.nBins;
            end
            
            sP = sqrt(nP);
            
            if round(sP) == sP
                rP = sP;
                cP = sP;
            else
                LP = floor(sP);
                HP = ceil(sP);
                
                lp = nP/LP;
                hp = nP/HP;
                
                if round(lp) == lp
                    rP = LP;
                    cP = lp;
                elseif round(hp) == hp
                    rP = hp;
                    cP = HP;
                else
                    rP = LP;
                    cP = HP;
                end
            end
            
            for i = 1:obj.nBins
                colB = obj.binCol(i,:);
                
                subplot(rP,cP,i)
                
                if isBckgrnd
                    fill(X,Y,colB,'FaceAlpha',0.1,'LineStyle','none');
                    hold on
                end
                
                [cladeTaxa,nC,colG] = obj.binTaxa(i).totalCladeSplit;
                
                for j = 1:nC
                    x = cladeTaxa(j).scores(:,xPC);
                    y = cladeTaxa(j).scores(:,yPC);
                    nT = cladeTaxa(j).num;
                    
                    if nT > 2
                        ch = convhull(x,y);
                        fill(x(ch),y(ch),colG(j,:),'FaceAlpha', 0.5);
                        hold on
                    end
                    plot(x,y,'ko','MarkerSize',5,'MarkerFaceColor',colG(j,:));
                    hold on
                end
                
                xlim(xL);
                ylim(yL);
                
                title(obj.binName(i));
                
            end
        end
        
        function Draw(obj,xPC,yPC)
            
            taxa = obj.allTaxa;
            
            minX = min(taxa.scores(:,xPC));
            minY = min(taxa.scores(:,yPC));
            maxX = max(taxa.scores(:,xPC));
            maxY = max(taxa.scores(:,yPC));
            
            xB = (maxX - minX) / 20;
            yB = (maxY - minY) / 20;
            
            xL = [(minX - xB) (maxX + xB)];
            yL = [(minY - yB) (maxY + yB)];
            
            X = [xL(1) xL(1) xL(2) xL(2)];
            Y = [yL(1) yL(2) yL(2) yL(1)];
            
            time = zeros(1,obj.nBins);
            t = ones(taxa.num,obj.nBins);
            
            dT = (obj.bins(end) - obj.bins(1)) / obj.nBins;
            
            dT = dT * 0.025;
            
            for i = 1:obj.nBins
                x = obj.binTaxa(i).scores(:,xPC);
                y = obj.binTaxa(i).scores(:,yPC);
                
                col = obj.binCol(i,:);
                nT = length(x);
                
                time(i) = (obj.bins(i) + obj.bins(i+1)) / 2;
                t(:,i) = t(:,i) * time(i) * -1;
                fill3(X,Y,(t(1:4,i))+dT,col,'FaceAlpha',0.25,'LineWidth',1.1);
                hold on
                ch = convhull(x,y);
                fill3(x(ch),y(ch),t(ch,i),col);
                plot3(x,y,t(1:nT,i)-dT,'ko','MarkerSize',3,'MarkerFaceColor','w');
            end
            
            text(xL(1) - xB, yL(2) + xB, obj.bins(1)*-1,num2str(obj.bins(1)));
            
            for i = 1:obj.nBins
                text(xL(1) - 2*xB, yL(2) + 2*xB, time(i)*-1,num2str(time(i)),'HorizontalAlignment','right');
                text(xL(1) - xB, yL(2) + xB, obj.bins(i + 1)*-1,num2str(obj.bins(i + 1)));
            end
            xlim(xL);
            ylim(yL);
            g = gca;
            g.PlotBoxAspectRatio = [1 1 obj.nBins];
            view(-45,45);
            axis off
        end
        
        function [L1,yD,B] = SliceTransform(obj,xL,yL,r)
            ax = gca;
            pos = ax.Position;
            
            ctr = obj.binCtr;
            
            h = pos(4);
            w = pos(3);
            
            tR = obj.bins(end)-obj.bins(1);
            L = w/tR;
            
            B = (obj.bins - obj.bins(1)).*L;
            ctr = (ctr - obj.bins(1)).*L;
            
            xR = xL(2) - xL(1);
            yR = yL(2) - yL(1);
            
            if nargin < 4
                r = xR/yR;
            end
            
            a = (r^2)-1;
            b = 2*h;
            L1 = ctr(1) - B(1);
            c = (L1^2) + h^2;
            
            if a == 0
                yD = c/b;
            else
                yD = abs(((-1*b) + sqrt((b^2) - (4*a*c)))/(2*a));
            end
            
            B = horzcat(0,ctr);
        end
        
        function DrawSlice(obj,bin,L1,yD,B)
            ax = gca;
            h = ax.Position(4);
            
            Y = [0 yD h h-yD];
            X = [B(bin+1)-L1 B(bin+1)-L1 B(bin+1) B(bin+1)];
            
            fill(X,Y,obj.binCol(bin,:),'FaceAlpha',0.1,'LineStyle','none');
            hold on
            plot([X(4) X(4) X(1) X(1) X(4)],[Y(3) Y(4) Y(1) Y(2) Y(3)],'k-','LineWidth',1.1);
        end
        
        function [axList] = DrawGrid(obj,ax)
            if nargin < 2
                ax = axes;
            else
                axes(ax);
            end
            
            pos = ax.Position;
            
            h = pos(4);
            w = pos(3);
            
            axis off
            barAx = axes;
            barAx.Position = [pos(1) pos(2)+h*0.9 w h*0.1];
            obj.binBar(true);
            tAx = axes;
            tAx.Position = [pos(1) pos(2) w h*0.9];
            
            P1 = tAx.Position(1);
            P2 = tAx.Position(2);
            
            W = tAx.Position(3);
            H = tAx.Position(4);
            w = W/(obj.nBins);
            h = 4*H/10;
            
            G = H/10;
            B1 = 0;
            B2 = h;
            B3 = h*2;
            
            Y1 = [B1 B2 B2 B1];
            Y2 = [B2 B3 B3 B2] + G;
            Y3 = [B3+G H H B3+G];
            
            tr = obj.bins(end) - obj.bins(1);
            
            
            for i = 1:obj.nBins
                L = (i-1)*w;
                R = i*w;
                
                TL = W*(obj.bins(i) - obj.bins(1))/tr;
                TR = W*(obj.bins(i+1) - obj.bins(1))/tr;
                
                X = [L L R R];
                XT = [L TL TR R];
                
                fill(X,Y1,obj.binCol(i,:),'LineStyle','none','FaceAlpha',0.25);
                hold on
                fill(X,Y2,obj.binCol(i,:),'LineStyle','none','FaceAlpha',0.25);
                fill(XT,Y3,obj.binCol(i,:),'LineWidth',0.5,'FaceAlpha',0.1);
                
                axList(i,1) = axes('Position',[P1+L P2 w h],'Color','none');
                axList(i,2) = axes('Position',[P1+L P2+h+G w h],'Color','none');
                axes(tAx);
                fill(X,[B2 B2+G B2+G B2],'w','LineWidth',1.1,'FaceAlpha',0);
            end
            fill([0 0 W W], [B1 B3+G B3+G B1],'w','LineWidth',1.1,'FaceAlpha',0);
            fill([0 0 W W], [B1 H H B1],'w','LineWidth',1.1,'FaceAlpha',0);
            
            
            uistack(tAx,'bottom')
            xlim([0 W]);
            ylim([0 H]);
            axis off
        end
        
        function mD = MorphRankGrid(obj,POL,ax)
            
            
            if nargin < 3
                ax = axes;
            else
                axes(ax);
            end
            
            f = gcf;
            gridAx = obj.DrawGrid(ax);
            
            x = obj.allTaxa.scores(:,1);
            y = obj.allTaxa.scores(:,2);
            
            
            minPX = 0;
            maxPX = 1;
            
            
            xpbrd = (maxPX-minPX)*0.1;
            
            minMX = min(x);
            maxMX = max(x);
            minMY = min(y);
            maxMY = max(y);
            
            xmbrd = (maxMX-minMX)*0.1;
            ymbrd = (maxMY-minMY)*0.1;
            
            xPL = [minPX-xpbrd maxPX+xpbrd];
            
            xML = [minMX-xmbrd maxMX+xmbrd];
            yML = [minMY-ymbrd maxMY+ymbrd];
%             
%             xPR = xPL(2) - xPL(1);
%             
%             xMR = xML(2) - xML(1);
%             yMR = yML(2) - yML(1);
%             
%             Z = (Z - xPL(1))./xPR;
            mD = zeros(obj.nBins,1);
            cD = zeros(obj.nBins,2);
            
            for i = 1:obj.nBins
%                 x0 = (i-1)*w;
                
                
                X = obj.binTaxa(i).scores(:,1);
                Y = obj.binTaxa(i).scores(:,2);
                
                z = POL.getZ(X,Y);
                for j = 1:length(z)
                   if z > 1
                       z = 1;
                   elseif z < 0
                       z = 0;
                   end
                end
                ksx = [0:0.01:1];
                                
                phat = mle(z,'distribution','normal');
                pd = makedist('normal','mu',phat(1),'sigma',phat(2));
                pd = truncate(pd,0,1);
                D = pdf(pd,ksx);
                maxPY = max(D);
                minPY = min(D);
                ypbrd = (maxPY-minPY)*0.1;
                yPL = [minPY-ypbrd maxPY+ypbrd];
                disp(yPL);
                yR = yPL(2) - yPL(1);
                mD(i) = pd.mean;
                
                ch = convhull(X,Y);
                
                axes(gridAx(i,1));
                D = (D - minPY)./(maxPY-minPY);
                plot(ksx,D,'k-');
                hold on
                plot([mD(i) mD(i)],[0 1],'r-','LineWidth',1.1);
                xlim(xPL);
                ylim((yPL - minPY)./(maxPY-minPY));
                gridAx(i,1).Color = 'none';
                xT = [ceil(minPX*10)/10 floor(maxPX*10)/10];
                yT = [0 1];
                if i == 1
                    set(gridAx(i,1),'xtick',[]);
                    xlabel('Rank');
                    set(gridAx(i,1),'ytick',[]);
                    ylabel('Density');
                else
                    set(gridAx(i,1),'xtick',[]);
                    xlabel('Rank');
                    set(gridAx(i,1),'ytick',[]);
                end
                
                axes(gridAx(i,2));
                fill(X(ch),Y(ch),obj.binCol(i,:));
                hold on
                plot(X,Y,'ko','MarkerSize',1);
                xlim(xML);
                ylim(yML);
                gridAx(i,2).Color = 'none';
                
                xT = [ceil(minMX*10)/10 floor(maxMX*10)/10];
                yT = [ceil(minMY*10)/10 floor(maxMY*10)/10];
                
                if i == 1
                    set(gridAx(i,2),'xtick',[]);
                    xlabel('PC1');
                    set(gridAx(i,2),'ytick',[]);
                    ylabel('PC2');
                else
                    set(gridAx(i,2),'xtick',[]);
                    xlabel('PC1');
                    set(gridAx(i,2),'ytick',[]);
                end

            end
        end
        
        function pt = Transform(obj,x,y,bin,xL,yL,r)
            if nargin < 7
                [L1,yD,B] = obj.SliceTransform(xL,yL);
            else
                [L1,yD,B] = obj.SliceTransform(xL,yL,r);
            end
            
            x = x - xL(1);
            y = y - yL(1);
            
            xR = xL(2) - xL(1);
            yR = yL(2) - yL(1);
            
            x = x./xR;
            y = y./yR;
            
            ax = gca;
            h = ax.Position(4);
            
            X = B(bin+1)-L1;
            xV = [L1 h-yD];

            pt = zeros(length(x),2);

            for i = 1:length(x)
                pt(i,:) = xV*x(i);
                pt(i,1) = pt(i,1) + X;
                pt(i,2) = pt(i,2) + y(i)*yD;
            end
        end
        
        function DrawParetoSlice(obj,x,y,P)
            P.PlotTransform(x,y)
        end
        
        function PlotPoints(obj,x,y,col)
            plot(x,y,'ko','MarkerFaceColor',col,'MarkerEdgeColor',col,'MarkerSize',2);
        end
        
        function DrawMorphSlice(obj,x,y,Arg)
            % Arg  = [pcx, pcy, bin];
            
            pcx = Arg(1);
            pcy = Arg(2);
            bin = Arg(3);
            
            oX = obj.binTaxa(bin).scores(:,pcx);
            oY = obj.binTaxa(bin).scores(:,pcy);
            
            ch = convhull(oX,oY);
            
            fill(x(ch),y(ch),obj.binCol(bin,:),'FaceAlpha',0.8);
            plot(x,y,'ko','MarkerSize',1,'MarkerFaceColor','k');
        end
        
        function DrawPareto(obj,LS1,LS2,opt1,opt2,xPC,yPC,ax)
            f = figure;
            f.Units = 'pixels';
            f.Position(4) = f.Position(3);
            f.Units = 'normalized';
            f.Position(2) = 0.2;
            P = pareto(LS1.zvec,LS2.zvec,opt1,opt2);
            
            minX = min(P.v1);
            maxX = max(P.v1);
            minY = min(P.v2);
            maxY = max(P.v2);
            
            xbrd = (maxX-minX)*0.1;
            ybrd = (maxY-minY)*0.1;
            
            xL = [minX-xbrd maxX+xbrd];
            yL = [minY-ybrd maxY+ybrd];
            
            if nargin < 8
                ax = axes;
            else
                axes(ax);
            end
            
            pos = ax.Position;
            
            h = pos(4);
            w = pos(3);
            
            axis off
            barAx = axes;
            barAx.Position = [pos(1) pos(2)+h*0.95 w h*0.05];
            obj.binBar(false);
            tAx = axes;
            tAx.Position = [pos(1) pos(2) w h*0.95];
            
            
            for i = 1:obj.nBins
                x = obj.binTaxa(i).scores(:,xPC);
                y = obj.binTaxa(i).scores(:,yPC);
                v1 = LS1.getZ(x,y);
                v2 = LS2.getZ(x,y);
                [L1,yD,B] = obj.SliceTransform(xL,yL,1);
                Ppts = obj.Transform(P.v1,P.v2,i,xL,yL,1);
                TPpts = obj.Transform(v1,v2,i,xL,yL,1);
                obj.DrawSlice(i,L1,yD,B);
                hold on
                obj.DrawParetoSlice(Ppts(:,1),Ppts(:,2),P);
                obj.PlotPoints(TPpts(:,1),TPpts(:,2),[1 0 0]);
            end
            
            xlim([0 tAx.Position(3)]);
            ylim([0 tAx.Position(4)]);
            axis off
        end
        
        function DrawMorph(obj,xPC,yPC,ax)
            taxa = obj.allTaxa;
            minX = min(taxa.scores(:,xPC));
            maxX = max(taxa.scores(:,xPC));
            minY = min(taxa.scores(:,yPC));
            maxY = max(taxa.scores(:,yPC));
            
            brd = (maxY-minY)*0.1;
            
            xL = [minX-brd maxX+brd];
            yL = [minY-brd maxY+brd];
            
            if nargin < 4
                ax = axes;
            else
                axes(ax);
            end
            
            pos = ax.Position;
            
            h = pos(4);
            w = pos(3);
            
            axis off
            barAx = axes;
            barAx.Position = [pos(1) pos(2)+h*0.95 w h*0.05];
            obj.binBar(false);
            tAx = axes;
            tAx.Position = [pos(1) pos(2) w h*0.95];
            
            for i = 1:obj.nBins
                x = obj.binTaxa(i).scores(:,xPC);
                y = obj.binTaxa(i).scores(:,yPC);
                plotArg = [xPC yPC i];
                [L1,yD,B] = obj.SliceTransform(xL,yL);
                pts = obj.Transform(x,y,i,xL,yL);
                obj.DrawSlice(i,L1,yD,B);
                hold on
                obj.DrawMorphSlice(pts(:,1),pts(:,2),plotArg);
            end
            
            xlim([0 tAx.Position(3)]);
            ylim([0 tAx.Position(4)]);
            axis off
        end
        
        function ctr = binCtr(obj)
            ctr = zeros(1,obj.nBins);
            for i = 1:obj.nBins
                ctr(i) = (obj.bins(i) + obj.bins(i+1))/2;
            end
        end
        
        function plotTimeBand(obj,yL,col)
            if nargin < 3
                col = [0.9 0.9 0.9];
            end
            if nargin < 2
                yL = [-1 1];
            end
            for i = 1:obj.nBins
                X = [obj.bins(i) obj.bins(i) obj.bins(i+1) obj.bins(i+1)];
                Y = [yL(1) yL(2) yL(2) yL(1)];
                if mod(i,2)
                    fill(X,Y,col,'LineStyle','none');
                    hold on
                else
                    fill(X,Y,'w','LineStyle','none');
                    hold on
                end
            end
            xL = [obj.bins(1) obj.bins(end)];
            fill([xL(1) xL(1) xL(2) xL(2)],[yL(1) yL(2) yL(2) yL(1)],'w','FaceAlpha',0,'LineWidth',1.1);
            xlim(xL);
            ylim(yL);
        end
        
        %% DISPARITY
        
        function [m,ci] = kslim(obj,X,lim)
            [f,xi] = ksdensity(X);
            pd = fitdist(X,'Kernel');
            m = mean(pd);
            int = zeros(length(f),1);
            for i = 1:length(f)-1
                int(i+1) = (f(i+1) + f(i))*(xi(i+1) - xi(i))/2 + int(i);
            end
            int = int./int(end);
            for i = 1:length(f)-1
                if int(i) <= lim(1) && int(i+1) > lim(1)
                    ci(1) = xi(i);
                end
                if int(i) <= lim(2) && int(i+1) > lim(2)
                    ci(2) = xi(i);
                end
            end
        end
        
        function [D] = disparityCurve(obj,funct,varName,nSamp)
            UD = zeros(obj.nBins,1);
            MD = zeros(obj.nBins,1);
            LD = zeros(obj.nBins,1);
            d = zeros(nSamp,obj.nBins);
            for i = 1:obj.nBins
                X = obj.binTaxa(i).bootstrap(varName,nSamp);
                for j = 1:nSamp
                    d(j,i) = funct(X(:,:,j));
                end
                [m,ci] = obj.kslim(d(:,i),[0.05 0.95]);
                MD(i) = m;
%                 s = std(d(:,i));
%                 UD(i) = MD(i) + 2*s;
%                 LD(i) = MD(i) - 2*s;
                UD(i) = ci(2);
                LD(i) = ci(1);
            end
            p = zeros(1,obj.nBins-1);
            for i = 1:obj.nBins-1
                [~,p(i)] = ttest(d(:,i),d(:,i+1));
            end
            H = obj.pHB(p);
            disp(H);
            D = horzcat(LD,MD,UD);
        end
        
        function plotCurve(obj,D,isBinCol,col)
            if nargin < 4
                col = [0.1 0.1 0.1];
            end
            ctr = obj.binCtr;
            xBnd = horzcat(ctr,fliplr(ctr));
            yBnd = vertcat(D(:,3),flipud(D(:,1)));
            if isBinCol
                brd = (max(D(:,3)) - min(D(:,1)))*0.02;
                yD = vertcat(D(:,2)-brd,flipud(D(:,2)+brd));
                vertexCol = vertcat(obj.binCol,flipud(obj.binCol));
                fill(xBnd,yBnd,col,'FaceVertexCData',vertexCol,'FaceColor','interp','EdgeColor',[0 0 0],'FaceAlpha',1);
                hold on
                plot(ctr,D(:,2),'ko-','MarkerFaceColor','w');
            else
                fill(xBnd,yBnd,col,'FaceAlpha',0.25,'EdgeColor',col,'EdgeAlpha',0.5);
                hold on
                plot(ctr,D(:,2),'k-','Color',col,'LineWidth',2);
            end
        end
        
        function lim = limCalc(obj,D)
            Min = min(min(D));
            Max = max(max(D));
            rng = Max - Min;
            brd = rng*0.1;
            lim = [Min-brd Max+brd];
        end
        
        function [X,lim] = scale(obj,D)
            lim = obj.limCalc(D);
            r = lim(2) - lim(1);
            [R,C] = size(D);
            X = zeros(R,C);
            for i = 1:R
                for j = 1:C
                    X(i,j) = (D(i,j) - lim(1))/r;
                end
            end
            
            X = X-0.5;
        end
        
        function binBar(obj,showNames)
            if nargin < 2
                showNames = true;
            end
            t = text(0.5, 0.5, 'Q','FontSize',8);
            E = t.Extent;
            Qw = E(3);
            Y = [0 1 1 0];
            ctr = obj.binCtr;
            TR = obj.bins(end) - obj.bins(1);
            C = (ctr - obj.bins(1))./TR;
            B = (obj.bins - obj.bins(1))./TR;
            for i = 1:obj.nBins
                X = [B(i) B(i) B(i+1) B(i+1)];
                fill(X,Y,obj.binCol(i,:),'LineWidth',1.1);
                hold on
                if showNames
                    n = length(obj.binName{i});
                    if n*Qw >= (X(3) - X(1))
                        if 4*Qw >= (X(3) - X(1))
                            S = obj.binName{i};
                            s = [S(1),'.'];
                            text(C(i),0.5,s,'HorizontalAlignment','center','FontSize',8);
                        else
                            S = obj.binName{i};
                            s = [S(1:3),'.'];
                            text(C(i),0.5,s,'HorizontalAlignment','center','FontSize',8);
                        end
                    else
                        text(C(i),0.5,obj.binName{i},'HorizontalAlignment','center','FontSize',8);
                    end
                end
            end
            xlim([0 1]);
            ylim([0 1]);
            axis off
        end
        
        function dispPlot(obj,nSamp,POL)
            [V] = obj.disparityCurve(@soVar,'harm',nSamp);
            [R] = obj.disparityCurve(@soRge,'harm',nSamp);
            [D] = obj.disparityCurve(@meanDist,'harm',nSamp);
            [Vol] = obj.disparityCurve(@Volume,'scores',nSamp);
            
            
            
            figure
            obj.plotTimeBand([0 1]);
            ax = gca;
            P = ax.Position;
            ax.Position = [P(1) P(2)+P(4)*0.35 P(3) P(4)*0.65];
            
            colormap hsv
            cmap = colormap;
            [nC,~] = size(cmap);
            rC = round(nC/4);
            col = zeros(4,3);
            for i = 1:4
%                 col(i,:) = cmap(rC*(i-1) + 1,:);
                col(i,:) = [rand rand rand];
                col(i,:) = col(i,:)/norm(col(i,:));
                HSV = rgb2hsv(col(i,:));
                HSV(2) = 1;
                HSV(3) = 0.7;
                col(i,:) = hsv2rgb(HSV);
            end
            axis off
            
            ax = gca;
            xL = ax.XLim;
            pos = ax.Position;
            
            h = pos(4)/4;
            
            ax1 = axes;
            ax1.Position = [pos(1) pos(2) pos(3) h];
            obj.plotCurve(V,true);
            xlim(xL);
            ylim(obj.limCalc(V));
            ylabel('Sum of Variances');
            ax1.Box = 'off';
            ax1.XColor = 'none';
            ax1.Color = 'none';
            
            ax2 = axes;
            ax2.Position = [pos(1) pos(2)+h pos(3) h];
            obj.plotCurve(R,true);
            ax2.YAxisLocation = 'right';
            xlim(xL);
            ylim(obj.limCalc(R));
            ylabel('Sum of Ranges');
            ax2.Box = 'off';
            ax2.XColor = 'none';
            ax2.Color = 'none';
            
            ax3 = axes;
            ax3.Position = [pos(1) pos(2)+2*h pos(3) h];
            obj.plotCurve(D,true);
            xlim(xL);
            ylim(obj.limCalc(D));
            ylabel('Mean Pairwise Distance');
            ax3.Box = 'off';
            ax3.XColor = 'none';
            ax3.Color = 'none';
            
            tAx = axes;
            tAx.Position = [P(1) P(2) P(3) P(4)*0.35];
            mD = obj.MorphRankGrid(POL,tAx);
            
            ax4 = axes;
            ax4.Position = [pos(1) pos(2)+3*h pos(3) h];
            obj.plotCurve(Vol,true);
            ax4.YAxisLocation = 'right';
            xlim(xL);
            ylim(obj.limCalc(Vol));
            ylabel('PCo Area');
            ax4.Box = 'off';
            ax4.XColor = 'none';
            ax4.Color = 'none';

%             obj.DrawMorph(1,2,tAx);
%             ctr = obj.binCtr;
%             t = ctr(1)-obj.bins(1);
%             T = obj.bins(end)-obj.bins(1);
%             xS = t*P(3)/T;
%             barAx2 = axes;
%             barAx2.Position = [P(1)-xS P(2) P(3) P(4)/40];
%             obj.binBar(false);

            figure
            plot(V(:,2),mD,'ko');
            [rho,pval] = corr(V(:,2),mD);
            disp(rho);
            disp(pval);
            figure
            plot(R(:,2),mD,'ko');
            [rho,pval] = corr(R(:,2),mD);
            disp(rho);
            disp(pval);
            figure
            plot(D(:,2),mD,'ko');
            [rho,pval] = corr(D(:,2),mD);
            disp(rho);
            disp(pval);
            figure
            plot(Vol(:,2),mD,'ko');
            [rho,pval] = corr(Vol(:,2),mD);
            disp(rho);
            disp(pval);
            
            div = D./mD;
            figure
            plot(obj.binCtr,div,'k-');
            
            N = zeros(obj.nBins,1);
            for i = 1:obj.nBins
                N(i) = obj.binTaxa(i).num;
            end
            
            [rho,pval] = corr(N,mD);
            disp(rho);
            disp(pval);
        end
        
        function RateOfChange(obj,nSamp,POL,pcx,pcy)
            [V] = obj.disparityCurve(@Volume,'scores',nSamp);
            ctr = obj.binCtr;
            
            dV = zeros(length(V)-1,1);
            for i = 1:length(V)-1
                M = (V(i+1) - V(i))/(ctr(i+1) - ctr(i));
                dV(i) = M/V(i);
            end
            
            figure
            plot(ctr(1:length(dV)),dV,'k-');
            
            figure
            R = zeros(obj.nBins);
            for i = 1:obj.nBins
                scores = obj.binTaxa(i).scores;
                x = scores(:,pcx);
                y = scores(:,pcy);
                ch = convhull(x,y);
                z = POL.getZ(x(ch),y(ch));
                m = mean(z);
                poly = polyshape(x(ch),y(ch));
                p = perimeter(poly);
                R(i) = p/m;
            end
            plot(ctr,R,'k-');
            
        end
        
        function N = TaxonBinSizes(obj)
            for i = 1:obj.nBins
                disp(obj.binName{i});
                disp(obj.binTaxa(i).num);
                N(i) = obj.binTaxa(i).num;
            end
        end
        
        function H = pHB(obj,p)
            [B,I] = sort(p);
            m = length(p);
            k = 0;
            for i = 1:m
                BM = 0.05/(m+1-i);
                if B(i) > BM
                    k = i;
                    break;
                end
            end
            H = zeros(1,m);
            if k == 0
                H = H+1;
            elseif k ~= 1
                H(I(1:k-1)) = 1;
            end
        end
        
        function TimePF(obj,pcx,pcy,MS,LS1,LS2,opt1,opt2)
            figure
            for i = 1:obj.nBins
                subplot(3,3,i);
                MS.taxaPF(LS1,LS2,opt1,opt2,obj.binTaxa(i).scores(:,pcx),obj.binTaxa(i).scores(:,pcy));
            end
        end

    end
end

