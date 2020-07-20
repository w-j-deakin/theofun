classdef Phylo
    
    properties
        taxa
        edges
        
    end
    
    methods
        function obj = Phylo(taxaSet,edgeFile,labelFile)
            obj.edges = importdata(edgeFile);
            obj.edges = obj.edges.data;
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
        
        function Plot(obj,minLength)
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
            nodes = PhyloNode.empty(0,1);
            
            [nR, ~] = size(obj.edges);
            
            for i = 1:nR+1
                nodes(i) = PhyloNode(i);
            end
            
            for i = 1:nR
                nodes(obj.edges(i,1)) = nodes(obj.edges(i,1)).AddChild(obj.edges(i,2));
                nodes(obj.edges(i,2)).parent = obj.edges(i,1);
            end
            
            anc = importdata(ancestorFile);
            anc = anc.data(1:end, 4:end);
            disp(anc);
            
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
        
    end
end

