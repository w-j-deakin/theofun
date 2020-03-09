classdef ele2D 
    
    properties
        nodes
    end
    
    methods
        function obj = ele2D(nodeIDX)
            nD = ndims(nodeIDX);
            nN = size(nodeIDX);
            if nD > 2
                error('multi-dimensional node index in element generation');
            elseif nN(1) > 1 && nN(2) > 1
                error('multi-dimensional node index in element generation');
            elseif nN(1) ~= 3 && nN(2) ~= 3
                error('2D element requires 3 vertices');
            else
                obj.nodes = nodeIDX;
            end
        end
        
         function tf = inPoints(obj,nP)
             maxN = max(obj.nodes);
             tf = maxN <= nP;
         end
         
         function tf = test_ele2D(obj)
             tf = true;
         end
         
    end
end

