classdef PhyloNode
    
    properties
        children
        parent
        index
        mindate
        maxdate
        y
    end
    
    methods
        function obj = PhyloNode(children,parent,index,mindate,maxdate,phylogeny)
            if nargin < 1
                obj.children = [];
                obj.parent = [];
                obj.index = [];
                obj.mindate = [];
                obj.maxdate = [];
                obj.y = [];
            else
                if isempty(children)
                    obj.children = [];
                else
                    obj.children = children;
                end
                obj.parent = parent;
                obj.index = index;
                if isempty(mindate)
                    obj.mindate = [];
                else
                    obj.mindate = mindate;
                end
                if isempty(maxdate)
                    obj.maxdate = [];
                else
                    obj.maxdate = maxdate;
                end
                obj.y = obj.UpdateY(phylogeny);
            end
        end
        
        function obj = AddChild(obj,child)
            if isempty(obj.children)
                obj.children = child;
            else
                add = true;
                for i = 1:length(obj.children)
                    if obj.children(i) == child
                        add = false;
                    end
                end
                if add
                    obj.children = horzcat(obj.children,child);
                end
            end
        end
       
        
        function objP = UpdateParent(obj,phylogeny,minL)
            objP = phylogeny.nodes(obj.parent);
            objP = objP.AddChild(obj.index);
            objP = objP.UpdateY(phylogeny);
            if isempty(objP.mindate)
               objP.mindate = obj.mindate;
               objP.maxdate = obj.mindate;
            elseif obj.mindate < objP.mindate
               objP.mindate = obj.mindate;
               objP.maxdate = obj.mindate;
            end
            if obj.maxdate-minL < objP.mindate
               objP.mindate = obj.mindate - minL;
               objP.maxdate = obj.mindate - minL;
            end
        end
        
        function Y = GetChildrenY(obj,phylogeny)
            Y = [];
            for i = 1:length(obj.children)
                Y(i) = phylogeny.nodes(obj.children(i)).y;
            end
        end
        
        function C = GetChildrenClade(obj,phylogeny)
            for i = 1:length(obj.children)
                if obj.children(i) <= length(phylogeny.taxa)
                    C(i) = phylogeny.taxa(obj.children(i)).clade;
                else
                    c = phylogeny.nodes(obj.children(i)).GetChildrenClade(phylogeny);
                    C(i) = c(1);
                    for j = 1:length(c)-1
                        if c(j) ~= c(j+1) || c(j) == 0 || c(j+1) == 0
                            C(i) = 0;
                        end
                    end
                end
            end
        end
        
        function D = GetChildrenDates(obj,phylogeny)
            for i = 1:length(obj.children)
                D(i) = phylogeny.nodes(obj.children(i)).mindate;
            end
        end
        
        function obj = UpdateY(obj,phylogeny)
            if obj.index > length(phylogeny.taxa)
                Y = obj.GetChildrenY(phylogeny);
                obj.y = mean(Y);
            else
                obj.y = -1*obj.index;
            end
        end
        
        function N = NumParents(obj,phylogeny)
            if isempty(obj.parent)
                N = 0;
            else
                N = phylogeny.nodes(obj.parent).NumParents(phylogeny) + 1;
            end
        end
        
        function tf = ShareParent(obj1,obj2)
            tf = obj1.parent == obj2.parent;
        end
        
        
        function Draw(obj,phylogeny,showRange,showClade)
            p = phylogeny.nodes(obj.parent);
            if showClade
                C = phylogeny.NodeColor(obj.index);
            else
                C = [0 0 0];
            end
            plot([obj.maxdate p.maxdate],[obj.y obj.y],'k-','Color',C);
            hold on
            if obj.mindate ~= obj.maxdate && showRange
                plot([obj.maxdate obj.mindate],[obj.y obj.y],'k-','Color',C,'LineWidth',2);
            end
            if ~isempty(obj.children)
                Y = obj.GetChildrenY(phylogeny);
                minY = min(Y);
                maxY = max(Y);
                
                plot([obj.maxdate obj.maxdate],[minY maxY],'k-','Color',C);
            end
        end
    end
end

