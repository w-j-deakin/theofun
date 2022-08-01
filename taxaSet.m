classdef taxaSet
    
    % Taxon EFA Dataset. Use 'taxaSetFromFile' to load data.
    
    properties
        num
        % Number of taxa
        taxa = taxonData.empty(0,1);
        % taxonData for each taxon
        nHarm
        % Optimal number of harmonics to be used in analyses
        scores
        % PCA scores
        coeff
        % PCA coefficients
        mu
        % Mean
        explained
        % Explained variance by each PC axis
        timeRange
        % [min time bin, max time bin]
        clades
        % [clade numbers]
        cNames
        % List of clade names
    end
    
    methods
        %% CONSTRUCTOR
        function obj = taxaSet(inTaxa,NH,inScr,inCo,inExp,inMu,cladeNames)
            obj.taxa = inTaxa;
            obj.nHarm = NH;
            obj.scores = inScr;
            obj.coeff = inCo;
            obj.mu = inMu;
            obj.explained = inExp;
            
            N = length(inTaxa);
            obj.num = N;
            
            clade = zeros(N,1);
            FAD = zeros(N,1);
            LAD = zeros(N,1);
            
            for i = 1:N
                clade(i) = inTaxa(i).clade;
                FAD(i) = inTaxa(i).FAD;
                LAD(i) = inTaxa(i).LAD;
            end
            
            obj.clades = unique(clade);
            obj.timeRange = [min(FAD) max(LAD)];
            
            if nargin < 7
                cladeNames = clades;
            end
            
            obj.cNames = cladeNames;
        end
        
        %% CATEGORISE TAXASET
        
        % return taxaset in a time bin
        function obj = timeSplit(inObj,time)
            
            if time < inObj.timeRange(1) || time > inObj.timeRange(2)
                obj = taxaSet.empty;
            else
                Coeff = inObj.coeff;
                nH = inObj.nHarm;
                m = inObj.mu;
                exp = inObj.explained;

                k = 0;
                N = inObj.num;
                for i = 1:N
                    if inObj.taxa(i).checkTime(time)
                        k = k+1;
                        Taxa(k) = inObj.taxa(i);
                        Scores(k,:) = inObj.scores(i,:);
                    end
                end

                obj = taxaSet(Taxa,nH,Scores,Coeff,exp,m,inObj.cNames);
            end
        end
        
        % test whether any taxa are ungrouped
        function tf = UnGrouped(obj)
            tf = false;
            for i = 1:obj.num
                if obj.taxa(i).clade == 0 || isnan(obj.taxa(i).clade) || isempty(obj.taxa(i).clade)
                    tf = true;
                    break;
                end
            end
        end
        
        % return taxaset in a clade
        function obj = cladeSplit(inObj,clade)
            
            inC = false;
            nC = length(inObj.clades);
            for i = 1:nC
                if inObj.clades(i) == clade
                    inC = true;
                    cName = inObj.cNames{i};
                end
            end
            
            if inC
                Coeff = inObj.coeff;
                nH = inObj.nHarm;
                exp = inObj.explained;
                m = inObj.mu;

                k = 0;
                N = inObj.num;
                for i = 1:N
                    if inObj.taxa(i).checkClade(clade)
                        k = k+1;
                        Taxa(k) = inObj.taxa(i);
                        Scores(k,:) = inObj.scores(i,:);
                    end
                end

                obj = taxaSet(Taxa,nH,Scores,Coeff,exp,m,cName);
            else
                obj = taxaSet.empty;
            end
        end
        
        % display clade sample sizes
        function DisplayCladeSize(obj)
            for i = 1:length(obj.clades)
                disp(obj.cNames{i});
                t = obj.cladeSplit(obj.clades(i));
                disp(t.num);
            end
        end
        
        % split taxaset based on all clades
        function [objs,nC,col] = totalCladeSplit(inObj)
            nC = length(inObj.clades);
            
            colormap hsv;
            cmap = colormap;
            
            nCol = length(cmap);
            
            rC = nCol/nC;
            
            col = zeros(nC,3);
            
            objs = taxaSet.empty(0,1);
            
            for i = 1:nC
                objs(i) = inObj.cladeSplit(inObj.clades(i));
                col(i,:) = cmap(round(1+(i-1)*rC),:);
            end
        end
        
        
        % bootstrap data
        function X = bootstrap(obj,varName,nSamp)
            if strcmp(varName,'scores')
                input = obj.scores(:,1:2);
            elseif strcmp(varName,'harm')
                input = zeros(obj.num,4*obj.nHarm-3);
                for i = 1:obj.num
                    input(i,:) = obj.taxa(i).harmRow(obj.nHarm);
                end
            else
                S = ['"',varName,'" is not a valid variable name'];
                error(S);
            end
            
            [nR,nC] = size(input);
            
            X = zeros(nR,nC,nSamp);
            for i = 1:nSamp
                r = randi(nR,nR,1);
                X(:,:,i) = input(r,:);
            end
        end
        
        %% PLOTS
        
        % test and plot morphosapces robustness to sample size
        function Robustness(obj,nSamp,nC)
            r = obj.num/20;
            line = zeros(20,nC);
            lim = zeros(20,nC);
            
            C = obj.coeff;
            
            for i = 1:nC
                C(:,i) = C(:,i)/norm(C(:,i));
            end
            
            for i = 1:20
                N(i) = round(r*i);
                diff = zeros(nSamp,nC);
                for j = 1:nSamp
                    p = randi(obj.num,N(i),1);
                    T = obj.taxa(p);
                    X = T(1).harmRow(obj.nHarm);
                    for k = 2:N(i)
                        X = vertcat(X,T(k).harmRow(obj.nHarm));
                    end
                    [c,~,~,~,~,m] = pca(X);
                    
                    [~,nc] = size(c);
                    for k = 1:nC
                        if k <= nc
                            c(:,k) = c(:,k)/norm(C(:,k));
                            diff(j,k) = AngleND(C(:,k),c(:,k));
                            if diff(j,k) > pi/2
                                diff(j,k) = pi - diff(j,k);
                            end
                        else
                            diff(j,k) = nan;
                        end
                    end
                    diffm(j) = norm(obj.mu - m);
                end
                line(i,:) = mean(diff);
                lim(i,:) = 2.576 * std(diff) / (sqrt(nSamp));
                linem(i) = mean(diffm);
                limm(i) = 2.576 * std(diffm) / (sqrt(nSamp));
            end
            figure
            errorbar(N,linem,limm);
            xlabel('Subsample Size');
            ylabel('Error in Mean');
            
            figure
            for i = 1:nC
                S = ['PC', num2str(i)];
                errorbar(N,line(:,i),lim(:,i),'DisplayName', S);
                hold on
            end
            ylim([0 pi/2]);
            legend
            xlabel('Subsample Size');
            ylabel('PC Orientation Difference (Radians)');
        end
        
        % scree plot of explained PCA variance
        function pcaScree(obj)

            n = length(obj.explained);
            x = 1:n;
            y = zeros(n,1);
            y(1) = obj.explained(1);
            
            for i = 2:n
                y(i) = y(i-1) + obj.explained(i);
            end
            
            plot(x,y,'ko-','MarkerFaceColor','w');
            xlabel('Number of PC axes');
            ylabel('Cumulative Explained Variance (%)');
        end
        
        % clade convex hull morphospace
        function Cmorph(obj,ax1PC,ax2PC,isFilled)
            
            if nargin < 4
                isFilled = true;
            end

            Sym = {'ko'; 'ks'; 'kd'; 'k^'; 'kv'};
            Sz = [5 6 6 5 5];
            [group, nC, col] = obj.totalCladeSplit;
            
            for i = 1:nC
                x = group(i).scores(:,ax1PC);
                y = group(i).scores(:,ax2PC);
                
                if length(x) > 2
                    ch = convhull(x,y);
                
                    if isFilled
                        fill(x(ch),y(ch),col(i,:),'FaceAlpha',0.2);
                    else
                        fill(x(ch),y(ch),col(i,:),'FaceAlpha',0);
                    end
                end
                hold on
            end
            
            ungrouped = obj.UnGrouped();
            
            if ungrouped
                pts(1) = plot(obj.scores(:,ax1PC),obj.scores(:,ax2PC),'ko','MarkerFaceColor','k','MarkerSize',3);
                legendNames = [{'Ungrouped'};obj.cNames];
            else
                legendNames = obj.cNames;
            end
            
            for i = 1:nC
                x = group(i).scores(:,ax1PC);
                y = group(i).scores(:,ax2PC);
                if ungrouped
                    idx = i + 1;
                else
                    idx = i;
                end
                pts(idx) = plot(x,y,Sym{mod(i-1,5)+1},'MarkerFaceColor',col(i,:),'MarkerSize',Sz(mod(i-1,5)+1));
                hold on
            end
            
            lgd = legend(pts,legendNames);
            legend('boxoff');
            lgd.FontSize = 10;
            lgd.Location = 'eastoutside';
            
            
            xL = ['PC', num2str(ax1PC), ' (', num2str(obj.explained(ax1PC)), '%)'];
            yL = ['PC', num2str(ax2PC), ' (', num2str(obj.explained(ax2PC)), '%)'];
            
            xlabel(xL);
            ylabel(yL);
            
        end
        
        % plot scores
        function Plot(obj,pcx,pcy,marker,col,isLabel)
            x = obj.scores(:,pcx);
            y = obj.scores(:,pcy);
            plot(x,y,marker,'MarkerFaceColor',col,'MarkerSize',3)
            
            if nargin < 6
                isLabel = true;
            end
            if isLabel
                xL = ['PC',num2str(pcx),' (',num2str(obj.explained(pcx)),'%)'];
                yL = ['PC',num2str(pcy),' (',num2str(obj.explained(pcy)),'%)'];
                xlabel(xL);
                ylabel(yL);
            end
        end
        
        % plot names
        function PlotNames(obj,pcx,pcy,size)
            ax = gca;
            rX = ax.XLim(2) - ax.XLim(1);
            rY = ax.YLim(2) - ax.YLim(1);
            xD = rX/obj.num;
            yD = rY/obj.num;
            
            for i = 1:obj.num
                hold on
                text(obj.scores(i,pcx)+xD,obj.scores(i,pcy)+yD,obj.taxa(i).name,'FontSize',size);
            end
        end
        
        %plot scores in 3D space
        function Plot3(obj,pcx,pcy,z,marker,col,isLabel)
            x = obj.scores(:,pcx);
            y = obj.scores(:,pcy);
            Z = ones(length(x));
            Z = Z*z;
            plot3(x,y,Z,marker,'MarkerFaceColor',col,'MarkerSize',3);
            if nargin < 7
                isLabel = true;
            end
            if isLabel
                xL = ['PC',num2str(pcx),' (',num2str(obj.explained(pcx)),'%)'];
                yL = ['PC',num2str(pcy),' (',num2str(obj.explained(pcy)),'%)'];
                xlabel(xL);
                ylabel(yL);
            end
        end
        
        % plot scores separated by clade
        function PlotClade(obj,pcx,pcy)
            [objs,nC,col] = obj.totalCladeSplit;
            Sym = {'ko'; 'ks'; 'kd'; 'k^'; 'kv'};
            for i = 1:nC
                objs(i).Plot(pcx,pcy,Sym{mod(i,5)+1},col(i,:),false);
                hold on
            end
        end
        
        % plot scores separated by clade in 3D space
        function PlotClade3(obj,pcx,pcy,z,col)
            if nargin < 5
                [objs,nC,col] = obj.totalCladeSplit;
            else
                [objs,nC,~] = obj.totalCladeSplit;
            end
            Sym = {'ko'; 'ks'; 'kd'; 'k^'; 'kv'};
            for i = 1:nC
                x = objs(i).scores(:,pcx);
                y = objs(i).scores(:,pcy);
                Z = ones(length(x),1);
                Z = Z*z;
                ch = convhull(x,y);
                fill3(x(ch),y(ch),Z(ch),col(i,:),'FaceAlpha',0.1,'LineWidth',2)
                hold on
                plot3(x,y,Z+1,Sym{i},'MarkerFaceColor',col(i,:));
            end
            
            xL = ['PC',num2str(pcx),' (',num2str(obj.explained(pcx)),'%)'];
            yL = ['PC',num2str(pcy),' (',num2str(obj.explained(pcy)),'%)'];
            xlabel(xL);
            ylabel(yL);
        end
        
        % generate and plot ShapeSpace
        function ss = theoMorph(obj,ax1PC,ax2PC,n,margin,showNames)
            
            if nargin < 6
                showNames = false;
            end
            
            x = obj.scores(:,ax1PC);
            y = obj.scores(:,ax2PC);
            maxX = max(x);
            maxY = max(y);
            minX = min(x);
            minY = min(y);
            
            xrange = maxX - minX;
            
            xtr = margin*xrange;
            
            xR = [(minX - xtr) (maxX + xtr)];
            yR = [(minY - xtr) (maxY + xtr)];
            
            ss = shapespace(obj,ax1PC,ax2PC,n,xR,yR);

            xR = ss.xR;
            yR = ss.yR;

            nx = ss.Nx;
            ny = ss.Ny;
            
            bX = (xR(2) - xR(1)) / (nx);
            bY = (yR(2) - yR(1)) / (ny);
            
            xL = [(xR(1) - bX) (xR(2) + bX)];
            yL = [(yR(1) - bY) (yR(2) + bY)];
            
            figure
            
            set(gcf,'defaultAxesDataAspectRatioMode','manual');
            set(gcf,'defaultAxesDataAspectRatio',[1 1 1]);
            
            ss.plotGrid();
            axis equal
            hold on
            obj.Cmorph(ax1PC,ax2PC,true);
            if showNames
                obj.PlotNames(ax1PC,ax2PC,7);
            end
            xlim(xL);
            ylim(yL);
            

            xT = [round(xR(1),2) 0 round(xR(2),2)];
            yT = [round(yR(1),2) 0 round(yR(2),2)];
            
            
            
            set(gca,'xtick',xT);
            set(gca,'xticklabel',{num2str(xT(1)),num2str(xT(2)),num2str(xT(3))});
            set(gca,'ytick',yT);
            set(gca,'yticklabel',{num2str(yT(1)),num2str(yT(2)),num2str(yT(3))});
            set(gcf,'defaultAxesDataAspectRatioMode','factory');
            set(gcf,'defaultAxesDataAspectRatio','factory');
            
        end
        
        % generate taxon Mesh List
        function ML = taxonML(obj,nE)
            for i = 1:obj.num
                harm = obj.taxa(i).harm(1:obj.nHarm,:);
                scr = obj.scores(i,:);
                shape(i) = theoShapeN(scr,harm);
            end
            ML = meshList(shape,nE);
        end
        
        % generate taxon Mesh List in specified PC axes
        function [ML,aML] = taxonAx(obj,nE,ax)
            for i = 1:obj.num
                apos = obj.scores(i,:);
                pos = zeros(1,length(apos));
                pos(ax) = obj.scores(i,ax);
                apos(ax) = 0;
                cf = obj.coeff;
                cf = cf';
                H = pos*cf + obj.mu;
                S(i) = theoShapeN(pos,H);
                aH = apos*cf + obj.mu;
                aS(i) = theoShapeN(apos,aH);
            end
            ML = meshList(S,nE);
            aML = meshList(aS,nE);
        end
        
        % generate random Mesh List
        function ML = randML(obj,n,nE,margin,isNorm)
            if nargin < 4
                margin = 0;
            end
            if nargin < 5
                isNorm = false;
            end
            [~,nA] = size(obj.scores);
            
            if isNorm
                s = std(obj.scores);
                s = s + s*margin;
            else
                r = zeros(1,nA);
                m = zeros(1,nA);
                b = zeros(1,nA);
                for i = 1:nA
                    m(i) = min(obj.scores(:,i));
                    r(i) = max(obj.scores(:,i)) - m(i);
                    b(i) = r(i)*margin;
                    m(i) = m(i) - b(i);
                    r(i) = r(i) + 2*b(i);
                end
            end
            
            S = theoShapeN.empty(0,1);
            
            figure
            set(gcf,'defaultAxesDataAspectRatioMode','manual');
            set(gcf,'defaultAxesDataAspectRatio',[1 1 1]);
            k = 0;
            while k < n
                if isNorm
                    pos = normrnd(0,s,1,nA);
                else
                    rn = rand(1,nA);
                    pos = (rn.*r) + m;
                end
                H = pos*obj.coeff' + obj.mu;
                shape = theoShapeN(pos,H);
                if ~shape.isIntersect
                    k = k + 1;
                    S(k) = shape;
                    disp(k);
                end
            end
            
            ML = meshList(S,nE);
        end
        
        % Count taxa
        function N = countTaxa(obj,pcx,pcy,x,y,r)
            bX = x-(r/2);
            tX = x+(r/2);
            
            bY = y-(r/2);
            tY = y+(r/2);
            
            N = 0;
            for i = 1:obj.num
                if obj.scores(i,pcx) >= bX && obj.scores(i,pcx) <= tX && obj.scores(i,pcy) >= bY && obj.scores(i,pcy) <= tY
                    N = N+1;
                end
            end
        end
        
        % For ancestral state reconstruction in R (geomorph)
        function GenerateHarmonicArray(obj,filename)
            A = zeros(obj.nHarm*2,2,length(obj.taxa));
            for i = 1:length(obj.taxa)
                for j = 1:obj.nHarm
                    A((j*2)-1,1,i) = obj.taxa(i).harm(j,1);
                    A((j*2)-1,2,i) = obj.taxa(i).harm(j,2);
                    A(j*2,1,i) = obj.taxa(i).harm(j,3);
                    A(j*2,2,i) = obj.taxa(i).harm(j,4);
                end
            end
            csvwrite(filename,A);
        end
        
    end
end

