classdef meshList
    
    % Vector of 2D mesh objects
    properties
        n
        mesh = mesh2D.empty(0,1)
        shapes = theoShapeN.empty(0,1)
    end
    
    methods
        %% CONSTRUCTOR
        function [obj,meshed] = meshList(shapes,nE)
            N = length(shapes);
            k = 0;
            p = zeros(N,1);
            for i = 1:N
                if ~shapes(i).checkTSN
                    k = k + 1;
                    obj.shapes(k) = shapes(i);
                    M = generateMesh2D(shapes(i).harm,nE);
                    if isempty(M)
                        obj.shapes(k) = [];
                        k = k - 1;
                        disp('!');
                    else
                        meshed(k) = i;
                        obj.mesh(k) = M;
                    end
                end
                LoadBar(i,N,0);
            end
            
            obj.n = k;
        end
        
        %% FUNCTIONS
        
        % add 2 meshLists together
        function obj = addML(obj1,obj2)
            obj = meshList.empty;
            obj(1).n = obj1.n + obj2.n;
            obj(1).mesh = vertcat(obj1.mesh,obj2.mesh);
            obj(1).shapes = vertcat(obj1.shapes,obj2.shapes);
            obj = obj(1);
        end
        
        % perform autp FEA on each mesh
        function vms = autoFEA(obj,t,E,v)
            if nargin < 2
                t = 0.01;
            end
            if nargin < 3
                E = 2*(10^9);
            end
            if nargin < 4
                v = 0.3;
            end
            
            nM = obj.n;
            vms = zeros(nM,1);
            p = zeros(nM,1);
            for i = 1:nM
                vms(i) = median(obj.mesh(i).FEA(t,E,v));
                if vms(i) == 0
                    vms(i) = nan;
                end
                p(i) = i*100/nM;
                if (p(i) >= 25 && p(i-1) < 25) || (p(i) >= 50 && p(i-1) < 50) || (p(i) >= 75 && p(i-1) < 75)
                    S = [num2str(p(i)),'%'];
                    disp(S);
                end
            end
        end
        
        % perform bootstrapped FEA on each mesh
        function vms = randFEA(obj,n,r,aR,t,E,v)
            if nargin < 3
                r = 0.1;
            end
            if nargin < 4
                aR = [-pi/4 pi/4];
            end
            if nargin < 5
                t = 0.01;
            end
            if nargin < 6
                E = 2*(10^9);
            end
            if nargin < 7
                v = 0.3;
            end
            
            nM = obj.n;
            vms = zeros(n,nM);
            for i = 1:nM
                vms(:,i) = median(obj.mesh(i).btstrpFEA(n,[i nM],r,aR,t,E,v));
            end
        end
        
        % perform auto RI on each mesh
        function [RI,v] = autoRI(obj)
            nM = obj.n;
            v = zeros(nM,1);
            RI = zeros(nM,1);
            p = zeros(nM,1);
            
            for i = 1:nM
                [v(i),RI(i)] = obj.mesh(i).rotV;
                if v(i) == 0
                    v(i) = nan;
                end
                if RI(i) == 0
                    RI(i) = nan;
                end
                p(i) = i*100/nM;
                if (p(i) >= 25 && p(i-1) < 25) || (p(i) >= 50 && p(i-1) < 50) || (p(i) >= 75 && p(i-1) < 75)
                    S = [num2str(p(i)),'%'];
                    disp(S);
                end
            end
        end
        
        % perform bootstrapped RI on each mesh
        function [RI,v] = randRI(obj,n,r)
            
            if nargin < 3
                r = 0.1;
            end
            
            nM = obj.n;
            v = zeros(n,nM);
            RI = zeros(n,nM);
            
            for i = 1:nM
                [v(:,i),RI(:,i)] = obj.mesh(i).btstrpRI(n,[i nM],r);
            end
        end
        
    end
end

