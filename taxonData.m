classdef taxonData
    
    % Data on an individual taxon
    
    properties
        harm
        pwr
        name
        clade
        FAD
        LAD
    end
    
    methods
        %% CONSTRUCTOR
        function obj = taxonData(inName,inFAD,inLAD,inClade,nLdmk)
            obj.name = inName;
            obj.FAD = inFAD;
            obj.LAD = inLAD;
            obj.clade = inClade;
            
            eName = [char(inName), '.tif'];
            ldmk = autoLdmk(eName,nLdmk);
            
            if isempty(ldmk)
                obj.harm = [];
                obj.pwr = [];
            else
                [obj.harm,~,~,Pwr] = efa1(nLdmk,ldmk,40,10);
                obj.pwr = Pwr./Pwr(39);
            end
        end
        
        %% DATA
        
        % vectorise harmonic data
        function HR = harmRow(obj,nHarm)
            HR = obj.harm(1,4);
            for i = 2:nHarm
                newCols = [obj.harm(i,1) obj.harm(i,2) obj.harm(i,3) obj.harm(i,4)];
                HR = horzcat(HR,newCols);
            end
        end
        
        % check if taxon is within a stage
        function inTime = checkTime(obj,time)
            inTime = obj.FAD <= time && obj.LAD >= time;
        end
        
        % check if taxon is in a clade
        function inClade = checkClade(obj,clade)
            inClade = obj.clade == clade;
        end
        
        %% PLOTS
        
        % draw taxon
        function Draw(obj,col,resolution)
            [x,y] = efaDraw(obj.harm,resolution);
            fill(x,y,col);
            axis equal
        end
    end
end

