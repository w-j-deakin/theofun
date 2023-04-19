classdef Phylo
    
    properties
        taxa
        edges
        
    end
    
    methods
        function obj = Phylo(taxaSet,edgeFile,labelFile)
            obj.edges = importdata(edgeFile);
            
            tips = importdata(labelFile);
            obj.taxa = taxonData.empty(0,1);
            
            if isnumeric(tips)
                for i = 1:length(tips)
                    obj.taxa(i) = taxaSet.taxa(tips(i));
                end
            else
                for i = 1:length(tips)
                    found = false;
                    j = 1;
                    while ~found && j <= taxaSet.num
                        if strcmp(tips{i},taxaSet.taxa(j).name)
                            found = true;
                            obj.taxa(i) = taxaSet.taxa(j);
                        end
                        j = j + 1;
                    end
                    
                    if ~found
                        S = ['ERROR: could not find ', tips{i}, ' in taxaSet'];
                        error(S);
                    end
                end
            end
            
        end
        
        function nodes = CalculateNodes(obj)
            nodes = PhyloNode.empty(0,1);
            
            [nR, ~] = size(obj.edges);
            
            for i = 1:nR+1
                nodes(i) = PhyloNode(i);
                if i <= length(obj.taxa)
                    nodes(i).y = i;
                end
            end
            
            for i = 1:nR
                nodes(obj.edges(i,1)) = nodes(obj.edges(i,1)).AddChild(obj.edges(i,2));
                nodes(obj.edges(i,2)).parent = obj.edges(i,1);
            end
            
            for i = 1:nR+1
                nodes(i).y = nodes(i).CalculateY(nodes);
            end
        end
        
        function Plot(obj,minLength)
            nodes = obj.CalculateNodes();
            
            for i = 1:nR+1
                nodes(i).Plot(nodes,minLength);
                hold on
            end
            
            for i = 1:length(obj.taxa)
                text(minLength * -0.5,i,obj.taxa(i).name, 'FontSize', 7);
            end
            
            ax = gca;
            xL = ax.XLim;
            xL(2) = minLength * 3;
            ax.XLim = xL;
            
            axis off
        end
        
        function Phylomorphospace(obj,TDS,ancestorFile,pcx,pcy)
            nodes = obj.CalculateNodes();
            
            anc = importdata(ancestorFile);
            anc = anc.data(1:end, 4:end);
            
            [nr,nc] = size(anc);
            
            ancScr = zeros(nr,nc);
            taxScr = zeros(length(obj.taxa),nc);
            
            for i = 1:nr
                ancScr(i,:) = (anc(i,:)-TDS.mu) * TDS.coeff;
            end
            
            for i = 1:length(obj.taxa)
                V = obj.taxa(i).harmRow(TDS.nHarm);
                taxScr(i,:) = (V-TDS.mu) * TDS.coeff;
            end
            
            scr = vertcat(taxScr,ancScr);
            
            for i = 1:nR+1
                if ~isempty(nodes(i).parent)
                    plot([scr(i,pcx) scr(nodes(i).parent,pcx)], [scr(i,pcy) scr(nodes(i).parent,pcy)], 'k-');
                    hold on
                end
            end
            
            TDS.PlotClade(1,2);
            
        end
        
        function AncestorPhylo(obj,ancestorFile,minL,s)
            nodes = obj.CalculateNodes();
            
            anc = importdata(ancestorFile);
            anc = anc(1:end, 4:end);
            
            [nr,~] = size(anc);
            
            for i = 1:nr
                ancShapes(i,1) = theoShapeN(1,anc(i,:));
            end
            
            for i = 1:length(obj.taxa)
                taxaShapes(i,1) = theoShapeN(1,obj.taxa(i).harm);
            end
            
            shapes = vertcat(taxaShapes,ancShapes);
            
            
            
            figure
            
            for i = 1:length(nodes)
                nodes(i).Plot(nodes,minL);
                hold on
                x = nodes(i).GetMaxBranchLength(nodes) * minL * -1;
                [x,y] = shapes(i).draw(100,x,nodes(i).y,minL*s);
                fill(x,y,[0.75 0.75 0.75]);
            end
            
            for i = 1:length(obj.taxa)
                text(minL * -0.5,i,obj.taxa(i).name, 'FontSize', 7);
            end
            
            ax = gca;
            xL = ax.XLim;
            xL(2) = minL * 3;
            ax.XLim = xL;
            
            axis equal
            axis off
            
        end
        
    end
end

