classdef PhyloNode
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        index
        children
        parent
        y
    end
    
    methods
        function obj = PhyloNode(index, children, parent)
            obj.index = index;
            if nargin > 1
                obj.children = children;
                obj.parent = parent;
            else
                obj.children = [];
                obj.parent = [];
            end
            obj.y = [];
        end
        
        function obj = AddChild(obj,childIndex)
            if isempty(obj.children)
                obj.children = childIndex;
            else
                obj.children = horzcat(obj.children, childIndex);
            end
        end
        
        function Y = CalculateY(obj,nodes)
            if isempty(obj.children)
                Y = obj.y;
            else
                for i = 1:length(obj.children)
                    yC(i) = nodes(obj.children(i)).CalculateY(nodes);
                end
                Y = mean(yC);
            end
        end
        
        function N = GetMaxBranchLength(obj,nodes)
            if isempty(obj.children)
                N = 1;
            else
                for i = 1:length(obj.children)
                    n(i) = nodes(obj.children(i)).GetMaxBranchLength(nodes);
                end
                N = max(n) + 1;
            end
        end
        
        function Plot(obj,nodes,minLength)
            x = obj.GetMaxBranchLength(nodes) * minLength * -1;
            hold on
            if ~isempty(obj.children)
                xL = [ x x ];
                for i = 1:length(obj.children)
                    Y(i) = nodes(obj.children(i)).y;
                end
                yL = [ min(Y) max(Y) ];
                
                plot(xL, yL, 'k-');
            end
            if ~isempty(obj.parent)
                X = nodes(obj.parent).GetMaxBranchLength(nodes)*minLength*-1;
                xL = [ x X ];
                yL = [ obj.y obj.y ];
                
                plot(xL, yL, 'k-');
            end
        end
    end
end

