classdef Phylogeny
    
    properties
        nodes
        taxa
        binDates
    end
    
    methods
        %% CONSTRUCTOR
        function obj = Phylogeny(nodes,taxa,binDates)
            if nargin < 1
                obj.nodes = [];
            else
                obj.nodes = nodes;
            end
            if nargin < 2
                obj.taxa = [];
            else
                obj.taxa = taxa;
            end
            if nargin < 3
                obj.binDates = [];
            else
                obj.binDates = binDates;
            end
        end
        
        % build a tress from a newick file
        function obj = BuildTree(obj,treeFile,taxaSet,dates,minL)
            textdata = fileread(treeFile);
            findEnd = false;
            for i = 1:length(textdata)
                if i <= length(textdata)-11
                    if strcmp(textdata(i:i+11),'BEGIN TREES;')
                        begini = i+12;
                        findEnd = true;
                    end
                end
                if i <= length(textdata)-2
                    if strcmp(textdata(i:i+1),'EN') && findEnd
                        endi = i-2;
                        findEnd = false;
                    end
                end
            end
            
            tree = textdata(begini:endi);
            for i = 1:length(tree)
                if strcmp(tree(i),'=')
                    begini = i+2;
                end
            end
            tree = tree(begini:length(tree)-2);
            i = 1;
            k = 0;
            T = 0;
            while i <= length(tree)
                if strcmp(tree(i),'(') || strcmp(tree(i),')') || strcmp(tree(i),',')
                    i = i+1;
                else
                    for j = 1:3
                        if strcmp(tree(i+j),')') || strcmp(tree(i+j),',')
                            k = k+1;
                            T(k) = str2num(tree(i:i+j-1));
                            add = num2str(k);
                            tree = horzcat(tree(1:i-1),add,tree(i+j:end));
                            i = i+length(add);
                            break;
                        end
                    end
                end
            end
            
            obj.taxa = taxaSet.taxa(T);
            
            run = true;
            pos = 0;
            K = length(T);
            text = tree;
            while run
                pos = pos+1;
                if strcmp(text(pos),'(')
                    start = pos;
                elseif strcmp(text(pos),')')
                    bracket = text(start+1:pos-1);
                    for i = 1:length(bracket)
                        if strcmp(bracket(i),',')
                            bracket(i) = ' ';
                        end
                    end
                    nds = str2num(bracket);
                    K = K+1;
                    for i = 1:length(nds)
                        Nodes(nds(i)) = PhyloNode();
                        Nodes(nds(i)).index = nds(i);
                        Nodes(nds(i)).parent = K;
                        if nds(i) <= length(T)
                            Nodes(nds(i)).mindate = dates(obj.taxa(nds(i)).FAD);
                            Nodes(nds(i)).maxdate = dates(obj.taxa(nds(i)).LAD + 1);
                        end
                    end
                    add = num2str(K);
                    text = horzcat(text(1:start-1),add,text(pos+1:end));
                    pos = 0;
                end
                if ~strcmp(text(1),'(')
                    run = false;
                end
            end
            Nodes = horzcat(Nodes,PhyloNode([],[],length(Nodes)+1,[],[],obj));
            clc;
            obj.binDates = dates;
            obj.nodes = Nodes;
            
            for i = 1:length(T)
                obj.nodes(i) = obj.nodes(i).UpdateY(obj);
            end
            
            for i = 1:length(obj.nodes)-1
                obj.nodes(obj.nodes(i).parent) = obj.nodes(i).UpdateParent(obj,minL);
            end
        end
        
        % attempt a simple time scale (NOTE - DO NOT USE THIS FOR ANY
        % ANLYSIS. THIS IS FOR VISUALISATION ONLY).
        function obj = TimeTree(obj,minL)
            for i = 1:length(obj.taxa)
                FAD = obj.taxa(i).FAD;
                LAD = obj.taxa(i).LAD + 1;
                obj.nodes(i).mindate = obj.binDates(FAD);
                obj.nodes(i).maxdate = obj.binDates(LAD);
            end
            
            for i = length(obj.taxa):length(obj.nodes)
                D = obj.nodes(i).GetChildrenDates(obj);
                obj.nodes(i).mindate = min(D) - minL;
                obj.nodes(i).maxdate = obj.nodes(i).mindate;
            end
        end
        
        % rearrange the tree for more space conserving plots
        function obj = Rearrange(inObj,minL)
            obj = inObj;
            for i = -length(obj.nodes):-(length(obj.taxa)+1)
                i = i*-1;
                for j = 1:length(obj.nodes(i).children)
                    N(j) = obj.NumChildren(obj.nodes(i).children(j));
                end
                [~,idx] = sort(N);
                N = [];
                obj.nodes(i).children = obj.nodes(i).children(idx);
            end
            
            t = obj.GetFirstTaxa;
            
            
            for i = 1:length(obj.taxa)
                obj.nodes(t).index = i; 
                obj.nodes(t).y = -1*i;
                Nodes(i) = obj.nodes(t);
                Taxa(i) = obj.taxa(t);
                run = true;
                p = obj.nodes(obj.nodes(t).parent);
                if i < length(obj.taxa)
                    while run
                        if p.children(end) == t
                            t = obj.nodes(t).parent;
                            p = obj.nodes(obj.nodes(t).parent);
                        else
                            run = false;
                        end
                    end
                    for j = 1:length(p.children)-1
                        if p.children(j) == t
                            t = obj.GetFirstTaxa(p.children(j+1));
                            break;
                        end
                    end
                end
            end
            obj.nodes = horzcat(Nodes,obj.nodes(length(Nodes)+1:end));
            obj.taxa = Taxa;
            for i = 1:length(obj.nodes)
                obj.nodes(i).children = [];
            end
            
            for i = 1:length(obj.nodes)-1
                obj.nodes(obj.nodes(i).parent) = obj.nodes(i).UpdateParent(obj,minL);
            end
            
        end
        
        
        function obj = UpdateChildren(obj,node,dY)
            if isempty(obj.nodes(node).children)
                obj.nodes(node).y = obj.nodes(node).y + dY;
            else
                for i = 1:length(obj.nodes(node).children)
                    obj = obj.UpdateChildren(obj.nodes(node).children(i),dY);
                end
                obj.nodes(node).y = obj.nodes(node).y + dY;
            end
        end
        
        function N = NumChildren(obj,node)
            N = 0;
            if ~isempty(obj.nodes(node).children)
                for i = 1:length(obj.nodes(node).children)
                    N = N + obj.NumChildren(obj.nodes(node).children(i));
                end
                N = N + length(obj.nodes(node).children);
            end
        end
        
        function T = GetAllSubTaxa(obj,node)
            if isempty(obj.nodes(node).children)
                T = obj.nodes(node).index;
            else
                T = 0;
                for i = 1:length(obj.nodes(node).children)
                    T = horzcat(T,obj.GetAllSubTaxa(obj.nodes(node).children(i)));
                end
                T(1) = [];
            end
        end
        
        function t = GetFirstTaxa(obj,node)
            if nargin < 2
                node = length(obj.nodes);
            end
            if isempty(obj.nodes(node).children)
                t = node;
            else
                t = obj.GetFirstTaxa(obj.nodes(node).children(1));
            end
        end
        
        function [mindates,maxdates] = FindTimes(obj)
            for i = 1:length(obj.nodes)
               mindates(i) = obj.nodes(i).mindate;
               maxdates(i) = obj.nodes(i).maxdate;
            end
        end
        
        %% PLOTS
        
        function PlotTimeBand(obj,minL,yL,col)
            if nargin < 4
                col = [0.9 0.9 0.9];
            end
            if nargin < 3
                yL = [-(length(obj.taxa)+1) 0];
            end
            [mindates,maxdates] = FindTimes(obj);
            
            for i = 1:length(obj.binDates)-1
                X = [obj.binDates(i) obj.binDates(i) obj.binDates(i+1) obj.binDates(i+1)];
                Y = [yL(1) yL(2) yL(2) yL(1)];
                if mod(i,2)
                    fill(X,Y,col,'LineStyle','none');
                    hold on
                else
                    fill(X,Y,'w','LineStyle','none');
                    hold on
                end
            end
            xL = [min(mindates)-minL obj.binDates(end)+10*minL];
            fill([xL(1) xL(1) xL(2) xL(2)],[yL(1) yL(2) yL(2) yL(1)],'w','FaceAlpha',0,'LineWidth',1.1);
            xlim(xL);
            ylim(yL);
            yticks([]);
            xticks([]);
        end
        
        function C = NodeColor(obj,Node)
            L = length(obj.taxa);
            for i = 1:L
                c(i) = obj.taxa(i).clade;
            end
            numC = max(c);
            if Node <= L
                H = (obj.taxa(Node).clade-1)/numC;
                C = hsv2rgb([H 0.8 0.9]);
            else
                childC = obj.nodes(Node).GetChildrenClade(obj);
                colour = true;
                for i = 1:length(childC)-1
                    if childC(i) ~= childC(i+1) || childC(i) == 0 || childC(i+1) == 0
                        colour = false;
                        break;
                    end
                end
                if colour
                    H = (childC(1)-1)/numC;
                    C = hsv2rgb([H 0.8 0.9]);
                else
                    C = [0 0 0];
                end
            end
        end
        
        function tf = IsChildOf(obj,child,parent)
            if child == length(obj.nodes)
                tf = false;
            elseif obj.nodes(child).parent == parent
                tf = true;
            else
                tf = obj.IsChildOf(obj.nodes(child).parent,parent);
            end
        end
        
        function Plot(obj,minL,idx,showRange,showClade)
            if nargin < 4
                showRange = false;
            end
            if nargin < 5
                showClade = false;
            end
            
            figure
            obj.PlotTimeBand(minL);
            hold on
            for i = 1:length(obj.nodes)
                obj.nodes(i).Draw(obj,showRange,showClade);
            end
            for i = 1:length(obj.taxa)
                text(obj.nodes(i).maxdate+minL,obj.nodes(i).y,obj.taxa(i).name,'FontSize',7);
            end
            
%             for i = 1:length(obj.nodes)-length(obj.taxa)
%                 text(obj.nodes(idx(i)+length(obj.taxa)).maxdate+minL,obj.nodes(idx(i)+length(obj.taxa)).y,num2str(i+length(obj.taxa)),'FontSize',7);
%             end
            
        end
        
        function A = GenerateHarmonicArray(obj,nHarm)
            A = zeros(nHarm*2,2,length(obj.taxa));
            for i = 1:length(obj.taxa)
                for j = 1:nHarm
                    A((j*2)-1,1,i) = obj.taxa(i).harm(j,1);
                    A((j*2)-1,2,i) = obj.taxa(i).harm(j,2);
                    A(j*2,1,i) = obj.taxa(i).harm(j,3);
                    A(j*2,2,i) = obj.taxa(i).harm(j,4);
                end
            end
        end
        
        function E = GenerateEdgeArray(obj)
            E = zeros(length(obj.nodes)-1,2);
            for i = 1:length(obj.nodes)-1
                E(i,2) = i;
                E(i,1) = obj.nodes(i).parent;
            end
        end
        
        function P = FindCladisticPath(obj,index1,index2)
            run = true;
            k = 0;
            P1 = 0;
            i = index1;
            while run
                k = k+1;
                P1(k) = i;
                i = obj.nodes(i).parent;
                if isempty(i)
                    run = false;
                end
            end
            run = true;
            k = 0;
            P2 = 0;
            i = index2;
            while run
                k = k+1;
                P2(k) = i;
                i = obj.nodes(i).parent;
                if isempty(i)
                    run = false;
                end
            end
            
            P = 0;
            
            for i = 1:length(P1)
                for j = 1:length(P2)
                    if P1(i) == P2(j)
                        testP = horzcat(P1(1:i),fliplr(P2(1:j-1)));
                        if P == 0
                            P = testP;
                        elseif length(testP) < length(P)
                            P = testP;
                        end
                    end
                end
            end
            
        end
        
        function [A,idx] = ReOrderAncestors(obj,ancestors,edges)
            
            E = obj.GenerateEdgeArray();
            [N,~] = size(E);
            nT = length(obj.taxa);
            idx = zeros(N-nT+1,1);
            while ~all(idx)
                for i = 1:N
                    if idx(edges(i,1)-nT) == 0
                        if edges(i,2) <= nT
                            idx(edges(i,1)-nT) = obj.nodes(edges(i,2)).parent;
                        elseif idx(edges(i,2)-nT) ~= 0
                            idx(edges(i,1)-nT) = obj.nodes(idx(edges(i,2)-nT)).parent;
                        end  
                    end
                end
            end
            
            idx = idx-nT;
            [nr,nc] = size(ancestors);
            A = zeros(nr,nc);
            for i = 1:length(idx)
                A(idx(i),:) = ancestors(i,:);
            end
            
        end
        
        function [D,X,Y] = PatristicDistance(obj,index1,index2,ancestors)
            anc = ancestors(:,4:end);
            [~,n] = size(anc);
            nT = length(obj.taxa);
            
            X = obj.taxa(index1).harmRow(30);
            X = X(1:n);
            
            Y = obj.taxa(index2).harmRow(30);
            Y = Y(1:n);
            
            P = obj.FindCladisticPath(index1,index2);
            D = 0;
            
            for i = 1:length(P)-1
                if i == 1
                    V1 = X;
                    V2 = anc(P(i+1)-nT,:);
                elseif i == length(P)-1
                    V1 = anc(P(i)-nT,:);
                    V2 = Y;
                else
                    V1 = anc(P(i)-nT,:);
                    V2 = anc(P(i+1)-nT,:);
                end
                
                d = norm(V1-V2);
                D = D + d;
            end
            
        end
        
        function [rho,p] = PheneticVsPatristic(obj,nHarm,ancestors)
            nT = length(obj.taxa);
            PaD = zeros(nT*(nT-1)/2,1);
            k = 0;
            for i = 1:nT
                X(i,:) = obj.taxa(i).harmRow(nHarm);
            end
            PhD = pdist(X);
            for i = 1:nT-1
                for j = i+1:nT
                    k = k+1;
                    PaD(k) = obj.PatristicDistance(i,j,ancestors);
                end
            end
            
            
            figure
            plot(PhD',PaD,'ko');
            [rho,p] = corr(PhD',PaD);
            disp(rho);
            disp(p);
        end
        
        function PlotAnc(obj,ancestors,scale)
            obj.Plot(2);
            hold on
            nT = length(obj.taxa);
            for i = 1:nT
                H = obj.taxa(i).harm;
                x = obj.nodes(i).maxdate;
                y = obj.nodes(i).y;
                
                [ox,oy] = efaDraw(H,100);
                ox = ox*scale + x;
                oy = oy*scale + y;
                
                fill(ox,oy,[0 0 1]);
            end
            
            [nN,nH] = size(ancestors);
            nHarm = nH/4;
            
            for i = 1:nN
                X = ancestors(i,:);
                H = zeros(nHarm,4);
                for j = 1:nHarm
                    H(j,:) = X(4*j-3:4*j);
                end
                x = obj.nodes(i+nT).maxdate;
                y = obj.nodes(i+nT).y;
                
                [ox,oy] = efaDraw(H,100);
                ox = ox*scale + x;
                oy = oy*scale + y;
                fill(ox,oy,[1 0 0]);
            end
            axis equal
        end
        
       
    end
end

