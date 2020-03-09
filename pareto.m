classdef pareto
    
    % 2D Pareto Front.
    
    properties
        v1
        % vector of variable 1
        v2
        % vector of variable 2
        opt1
        % optimality of v1*
        opt2
        % optimality of v2*
        D
        % index list of pareto optimal values
        isempty
        % boolean to check that object has data (for Pareto Ranking)
        
        
        %  *If opt = 1 (true), v is optimal when maximised
        %   If opt = 0 (false), v is optimal when minimised
    end
    
    methods
        %% CONSTRUCTOR
        function obj = pareto(v1,v2,opt1,opt2)
            
            if nargin < 3
                opt1 = true;
            end
            if nargin < 4
                opt2 = true;
            end
            
            obj.v1 = v1;
            obj.v2 = v2;
            obj.opt1 = opt1;
            obj.opt2 = opt2;
            
            n = length(v1);
            
            if opt1
                t1 = v1;
            else
                t1 = v1*-1;
            end
            
            if opt2
                t2 = v2;
            else
                t2 = v2*-1;
            end
            
            k = 0;
            for i = 1:n
                if ~isnan(v1(i))
                    isDom = true;
                    run = true;
                    j = 1;
                    while run
                        if ~isnan(v2(j))
                            if (t1(j) > t1(i) && t2(j) >= t2(i)) || (t1(j) >= t1(i) && t2(j) > t2(i))
                                isDom = false;
                                run = false;
                            end
                        end
                        j = j+1;
                        if j > n
                            run = false;
                        end    
                    end
                    if isDom
                        k = k + 1;
                        D(k) = i;
                    end
                end
            end
            
            if k == 0
                obj.isempty = true;
                obj.D = 0;
            else
                obj.isempty = false;
                obj.D = D;
            end
            
        end
        
        %% FUNCTIONS
        
        % find pareto optimal subset
        function dom = FindDominant(obj,idx)
            if obj.opt1
                V1 = obj.v1;
            else
                V1 = obj.v1 * -1;
            end
            
            if obj.opt2
                V2 = obj.v2;
            else
                V2 = obj.v2 * -1;
            end
            
            T1 = V1(idx);
            T2 = V2(idx);
            k = 0;
            for i = 1:length(obj.v1)
                if V1(i) > T1 && V2(i) > T2
                    k = k+1;
                    dom(k) = i;
                end
            end
            if k == 0
                dom = [];
            end
        end
        
        % add two pareto front sets together
        function obj = add(obj1,obj2)
            % check both inputs are pareto class, with same optimality
            if obj1.opt1 == obj2.opt1 && obj1.opt2 == obj2.opt2
                V1 = vertcat(obj1.v1,obj2.v1);
                V2 = vertcat(obj1.v2,obj2.v2);
                obj = pareto(V1,V2,obj1.opt1,obj1.opt2);
            else
                error('check inputs are both of pareto class with same optimality');
            end
        end
        
        % find pareto rank of each solution
        function [R,p] = Rank(obj)
            
            V1 = obj.v1;
            V2 = obj.v2;
            o1 = obj.opt1;
            o2 = obj.opt2;
            d = obj.D;
            
            R = zeros(length(V1),1);
            
            R(d) = 1;
            
            V1(d) = nan;
            V2(d) = nan;
            p(1) = pareto(V1,V2,o1,o2);
            
            run = true;
            k = 1;
            
            while run
                d = p(k).D;
                k = k + 1;
                R(d) = k;
                V1(d) = nan;
                V2(d) = nan;
                
                p(k) = pareto(V1,V2,o1,o2);
                if p(k).isempty
                    run = false;
                end
            end
            
            mR = max(R);
            
            for i = 1:length(R)
                if R(i) == 0
                    R(i) = mR + 1;
                end
                if isnan(obj.v1(i))
                    R(i) = nan;
                end
            end
            
        end
        
        % find the rank ratio of each solution
        function R = RankRatio(obj)
            R = zeros(length(obj.v1),1); 
            ap = pareto(obj.v1,obj.v2,~(obj.opt1),~(obj.opt2));

            r1 = obj.Rank;
            r2 = ap.Rank;


            for i = 1:length(obj.v1)
                R(i) = (r2(i) - 1)/(r1(i) + r2(i) - 2);
            end
        end
        
        %% PLOTS
        
        function [x,y] = Line(obj)
            V1 = obj.v1;
            V2 = obj.v2;
            d = obj.D;
            [x,iX] = sort(V1(d));
            y = zeros(length(d),1);
            for i = 1:length(d)
                y(i) = V2(d(iX(i)));
            end
        end
        
        function Plot(obj,col,area,aP)
            
            [dX,dY] = obj.Line;
            if nargin < 3
                area = false;
            end
            if nargin < 4
                aP = pareto(obj.v1,obj.v2,~(obj.opt1),~(obj.opt2));
            end
            [sX,sY] = aP.Line;
            
            fX = vertcat(dX,flipud(sX));
            fY = vertcat(dY,flipud(sY));
            
            if nargin < 2
                fill(fX,fY,[0 0 0],'LineStyle','none','FaceAlpha',0.1);
                hold on
                plot(obj.v1,obj.v2,'k.','MarkerSize',1);
                plot(dX,dY,'k-','LineWidth',1.5);
            elseif area
                fill(fX,fY,col,'LineWidth',0.5,'EdgeColor',col,'FaceAlpha',0.25);
                hold on
                plot(dX,dY,'k-','LineWidth',1.5,'Color',col);
            else
                plot(obj.v1,obj.v2,'ko','MarkerSize',3,'MarkerFaceColor',col,'MarkerEdgeColor',col*0.5);
            end
%             set(gca,'Xscale','log');
        end
        
        function PlotTransform(obj,v1,v2,col)
            if nargin < 4
                col = [0.1 0.1 0.1];
            end
            p = pareto(obj.v1,obj.v2,obj.opt1,obj.opt2);
            ap = pareto(obj.v1,obj.v2,~obj.opt1,~obj.opt2);
            p.v1 = v1;
            ap.v1 = v1;
            p.v2 = v2;
            ap.v2 = v2;
            p.Plot(col,ap);
        end

        function overPlot(obj1,obj2)
        
            obj1.Plot;
            hold on
            [dX,dY] = obj2.Line;
            aP = pareto(obj2.v1,obj2.v2,~(obj2.opt1),~(obj2.opt2));
            [sX,sY] = aP.Line;
            
            fX = vertcat(dX,flipud(sX));
            fY = vertcat(dY,flipud(sY));
            fill(fX,fY,[0.9 0 0],'LineStyle','none','FaceAlpha', 0.2);
            hold on
            plot(obj2.v1,obj2.v2,'ro','MarkerSize',1);
            plot(dX,dY,'r-','LineWidth',3);
            
        end
        
        function RankPlot(obj)
            [~,p] = obj.Rank;
            figure
            [x,y] = obj.Line;
            plot(x,y,'k-','LineWidth',3);
            hold on
            
            for i = 1:length(p)-1
                [x,y] = p(i).Line;
                plot(x,y,'k:');
            end
        end

        function overRankPlot(obj1,obj2)
            obj1.RankPlot;
            hold on
            [~,p] = obj2.Rank;
            [x,y] = obj2.Line;
            plot(x,y,'r-','LineWidth', 3);
            hold on
            
            for i = 1:length(p)-1
                [x,y] = p(i).Line;
                plot(x,y,'r:');
            end
        end
        
        function [R,fR] = randRank(obj,nSamp,sampSize)
            
            N = length(obj.v1);
            n = sampSize;
            r = nan(N,nSamp);
            
            for i = 1:nSamp
                p = randperm(N,n);
                s1 = obj.v1(p);
                s2 = obj.v2(p);
                rP = pareto(s1,s2,obj.opt1,obj.opt2);
                r(p,i) = rP.Rank;
            end
            
            R = nanmean(r,2);
            k = 0;
            for i = 1:N
                if ~isnan(obj.v1(i)) && ~isnan(obj.v2(i)) && ~isnan(R(i))
                    k = k + 1;
                    V1(k) = obj.v1(i);
                    V2(k) = obj.v2(i);
                    nR(k) = R(i);
                end
            end
            
            fR = scatteredInterpolant(V1',V2',nR');
        end
        
    end
end

