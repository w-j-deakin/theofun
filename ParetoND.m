classdef ParetoND
    
    properties
        var
        opt
        D
        isempty
    end
    
    methods
        %% CONSTRUCTOR
        function obj = ParetoND(varin,optin)
            
            obj.var = varin;
            obj.opt = optin;
            
            [nT,nV] = size(varin);
            
            for i = 1:nV
                if ~optin(i)
                    varin(:,i) = varin(:,i) * -1;
                end
            end
            
            IDX = 0;
            for i = 1:nT
                if any(isnan(varin(i,:)))
                    isDom = false;
                else
                    isDom = true;

                    for j = 1:nT
                        if j ~= i
                            temp = zeros(1,nV);
                            for k = 1:nV
                                if varin(i,k) < varin(j,k)
                                    temp(k) = 1;
                                end
                            end
                            if sum(temp) == nV
                                isDom = false;
                            end
                        end
                    end
                end
                
                if isDom
                    IDX = IDX + 1;
                    D(IDX) = i;
                end
            end
            
            if IDX == 0
                obj.isempty = true;
                obj.D = 0;
            else
                obj.isempty = false;
                obj.D = D;
            end
            
        end
        
        %% FUNCTIONS
        
        % Goldberg rank
        function R = Rank(obj,i,inR)
            if obj.isempty
                R = inR;
            else
                if nargin == 1
                    [nT,~] = size(obj.var);

                    inR = zeros(nT,1);
                    for j = 1:nT
                        if any(isnan(obj.var(j,:)))
                            inR(j) = nan;
                        end
                    end
                    i = 1;
                end
                inR(obj.D) = i;
                
                newVar = obj.var;
                newVar(obj.D,:) = nan;
                i = i + 1;
                
                p = ParetoND(newVar,obj.opt);
                R = p.Rank(i,inR);
            end
        end
        
        % Rank Ratio
        function R = RankRatio(obj)
            RO = obj.Rank;
            pRev = ParetoND(obj.var, ~obj.opt);
            RS = pRev.Rank;
            
            R = zeros(length(RO),1);
            
            for i = 1:length(RO)
                if isnan(RO(i))
                    R(i) = nan;
                elseif RO(i) == 1
                    R(i) = 1;
                else
                    R(i) = (RS(i) - 1) / (RO(i) + RS(i) - 2);
                end
            end
        end
        
        
    end
end

