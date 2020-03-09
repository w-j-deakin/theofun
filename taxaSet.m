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
            
            rC = round(nCol/nC);
            
            col = zeros(nC,3);
            
            objs = taxaSet.empty(0,1);
            
            for i = 1:nC
                objs(i) = inObj.cladeSplit(inObj.clades(i));
                col(i,:) = cmap(1+(i-1)*rC,:);
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
            Sz = [3 3 3 3 3];
            [group, nC, col] = obj.totalCladeSplit;
            
            for i = 1:nC
                x = group(i).scores(:,ax1PC);
                y = group(i).scores(:,ax2PC);
                
                ch = convhull(x,y);
                
                if isFilled
                    fill(x(ch),y(ch),col(i,:),'FaceAlpha',0.2);
                else
                    fill(x(ch),y(ch),col(i,:),'FaceAlpha',0);
                end
                hold on
                plot(x(ch),y(ch),Sym{i},'MarkerFaceColor',col(i,:),'MarkerSize',Sz(i));
            end
            
            for i = 1:nC
                x = group(i).scores(:,ax1PC);
                y = group(i).scores(:,ax2PC);
%                 pts(i) = plot(x,y,Sym{i},'MarkerFaceColor',col(i,:),'MarkerEdgeColor',col(i,:)*0.8,'MarkerSize',Sz(i));
                hold on
            end
            
%             lgd = legend(pts,obj.cNames);
%             legend('boxoff');
%             lgd.FontSize = 10;
%             lgd.Location = 'eastoutside';
            
            
%             xL = ['PC', num2str(ax1PC), ' (', num2str(obj.explained(ax1PC)), '%)'];
%             yL = ['PC', num2str(ax2PC), ' (', num2str(obj.explained(ax2PC)), '%)'];
%             
%             xlabel(xL);
%             ylabel(yL);
            
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
                objs(i).Plot(pcx,pcy,Sym{i},col(i,:),false);
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
        function ss = theoMorph(obj,ax1PC,ax2PC,n,margin,flip)
            
            if nargin < 6
                flip = false;
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
            
            if flip
                ss = shapespace(obj,ax2PC,ax1PC,n,yR,xR);
                xR = ss.yR;
                yR = ss.xR;
            
                nx = ss.Ny;
                ny = ss.Nx;
            else
                ss = shapespace(obj,ax1PC,ax2PC,n,xR,yR);
            
                xR = ss.xR;
                yR = ss.yR;
            
                nx = ss.Nx;
                ny = ss.Ny;
            end
            
            bX = (xR(2) - xR(1)) / (nx);
            bY = (yR(2) - yR(1)) / (ny);
            
            xL = [(xR(1) - bX) (xR(2) + bX)];
            yL = [(yR(1) - bY) (yR(2) + bY)];
            
            figure
            
            set(gcf,'defaultAxesDataAspectRatioMode','manual');
            set(gcf,'defaultAxesDataAspectRatio',[1 1 1]);
            
            ss.plotGrid(flip);
            axis equal
            hold on
            obj.Cmorph(ax1PC,ax2PC,true);
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
        
        
        % PhyloMorphospace
        
        function PhyloMorph(obj,pcx,pcy,phylogeny,ancestors)
            x = obj.scores(:,pcx);
            y = obj.scores(:,pcy);
            
%             plot(x,y,'ko','MarkerFaceColor','w','MarkerEdgeColor',[0.25 0.25 0.25],'MarkerSize',3);
%             hold on
            
            for i = 1:length(phylogeny.taxa)
                for j = 1:length(obj.taxa)
                    if strcmp(phylogeny.taxa(i).name,obj.taxa(j).name)
                        X(i) = x(j);
                        Y(i) = y(j);
                        break;
                    end
                end
            end
            
            
            for i = length(phylogeny.taxa)+1:length(phylogeny.nodes)
                Anc = ancestors(i-length(phylogeny.taxa),4:end);
                v = (Anc - obj.mu) * obj.coeff;
                X(i) = v(pcx);
                Y(i) = v(pcy);
            end
            
            for i = 1:length(phylogeny.nodes)-1
                plot([X(i) X(phylogeny.nodes(i).parent)], [Y(i) Y(phylogeny.nodes(i).parent)],'k-','Color',[0.25 0.25 0.25]);
                hold on
            end
%             for i = 1:length(phylogeny.taxa)
%                 C = phylogeny.NodeColor(i);
%                 plot(X(i),Y(i),'ko','MarkerFaceColor',C,'MarkerEdgeColor',[0.25 0.25 0.25],'MarkerSize',3);
%                 hold on
% %                text(X(i)+0.002,Y(i)+0.002,phylogeny.taxa(i).name,'FontSize',8);
%             end
%             obj.PlotClade(pcx,pcy);
            obj.Cmorph(pcx,pcy);
        end
        
    end
end

