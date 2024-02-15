classdef Phylo
    
    properties
        taxa
        edges
        ancestors
        edgeLengths
        nodes
    end
    
    methods
        %% CONSTRUCTOR
        function obj = Phylo(taxaSet,edgeFile,labelFile)
            obj.edges = importdata(edgeFile);
            
            tips = importdata(labelFile);
            obj.taxa = taxonData.empty(0,1);
            obj.ancestors = [];
            obj.edgeLengths = [];
            

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

            obj.nodes = obj.CalculateNodes();
            
        end

        %% ADD DATA
        % Import ancestral EFA data
        function obj = ImportEFAAncestors(obj,ancestorFile)
            anc = importdata(ancestorFile);
            anc = anc(1:end, 4:end);
            obj.ancestors = anc;
        end
        
        % Set edge length data
        function obj = SetEdgeLengths(obj,edgeLengths)
            [nr,nc] = size(edgeLengths);
            [ne,~] = size(obj.edges);
            if nr == ne && nc == 1
                obj.edgeLengths = edgeLengths';
            elseif nc == ne && nr == 1
                obj.edgeLengths = edgeLengths;
            else
                error('edgeLengths must be a vector of length num edges')
            end
        end
        
        %% INTERNAL FUNCTIONALITY
        % set the nodes
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

        % get parent edge lengths of nodes
        function EL = NodeEL(obj)
            if isempty(obj.edgeLengths)
                error('Phylo object has no edge length data.');
            end
            [~,sortID] = sort(obj.edges(:,2));
            nE = length(sortID);
            E = obj.edges(sortID,:);
            EL = obj.edgeLengths(sortID);
            for i = 2:nE
                if E(i,2) - E(i-1,2) == 2
                    EL = horzcat(EL(1:i-1),0,EL(i:end));
                end
            end
        end

        % get taxon morphological matrix
        function MT = TaxonMatrix(obj,nV)
            nT = length(obj.taxa);
            if nargin < 2    
                [nH,~] = size(obj.taxa(1).harm);
            else
                nH = (nV + 3)/4;
            end
            MT = zeros(nT,4*nH-3);
            for i = 1:nT
                MT(i,:) = obj.taxa(i).harmRow(nH);
            end
        end

        % find ancestral lineage of node
        function L = Lineage(obj,index)
            p = obj.nodes(index).parent;
            if isempty(p)
                L = index;
            else
                L = horzcat(index,obj.Lineage(p));
            end
        end

        % find MRCA of two nodes
        function [mrca,id1,id2] = MRCA(obj,L1,L2)
            if length(L1) == 1
                L1 = obj.Lineage(L1);
            end
            if length(L2) == 1
                L2 = obj.Lineage(L2);
            end
            mrca = 0;
            for i = 1:length(L1)
                for j = 1:length(L2)
                    if mrca == 0 && L1(i) == L2(j)
                        mrca = L1(i);
                        id1 = i;
                        id2 = j;
                    end
                end
            end
        end

        % calculate the pairwise total evolutionary morphological distance
        % between a set of nodes within the tree
        function [D,obj] = EvoMorphDistance(obj,nodeIDs, ancestorFile)
            % number of nodes to compare
            N = length(nodeIDs);
            % number of pairwise comparisons
            ND = N * (N-1) / 2;

            % output D will be an ND by 1 vector. Use squareform to convert
            % to N-byN matrix
            D = zeros(ND,1);

            % load ancestors
            if nargin < 3 && isempty(obj.ancestors)
                error('Phylo object does not contan ancestral data. Input a filename to import data');
            elseif nargin > 2
                obj = obj.ImportEFAAncestors(ancestorFile);
            end
            % combined taxon ancestor morphological matrix
            [~,nV] = size(obj.ancestors);
            M = vertcat(obj.TaxonMatrix(nV),obj.ancestors);

            id = 0;

            % for each pairwise combo, find MRCA and distance to that point
            for i = 1:N-1
                L1 = obj.Lineage(nodeIDs(i));
                for j = i+1:N
                    id = id + 1;
                    L2 = obj.Lineage(nodeIDs(j));
                    [~,id1,id2] = obj.MRCA(L1,L2);
                    for m = 1:id1-1
                        D(id) = D(id) + norm(M(L1(m),:)-M(L1(m+1),:));
                    end
                    for n = 1:id2-1
                        D(id) = D(id) + norm(M(L2(n),:)-M(L2(n+1),:));
                    end
                end
            end
        end

        % calculate the pairwise total branch length distance
        % between a set of nodes within the tree
        function [D,obj] = BranchLengthDistance(obj,nodeIDs,edgeLengths)
            % number of nodes to compare
            N = length(nodesIDs);
            % number of pairwise comparisons
            ND = N * (N-1) / 2;

            % output D will be an ND by 1 vector. Use squareform to convert
            % to N-byN matrix
            D = zeros(ND,1);

            % load ancestors
            if nargin < 3 && isempty(obj.edgeLengths)
                error('Phylo object does not contan edge length data. Input a filename to import data');
            elseif nargin > 2
                obj = obj.SetEdgeLengths(edgeLengths);
            end
            
            id = 0;

            EL = obj.NodeEL();

            % for each pairwise combo, find MRCA and distance to that point
            for i = 1:N-1
                L1 = obj.Lineage(nodeIDs(i));
                for j = i+1:N
                    id = id + 1;
                    L2 = obj.Lineage(nodeIDs(j));
                    [~,id1,id2] = obj.MRCA(L1,L2);
                    for m = 1:id1-1
                        if id1 == 1
                            disp('!!!!');
                        end
                        D(id) = D(id) + EL(m);
                    end
                    for n = 1:id2-1
                        if id2 == 1
                            disp('!!!!');
                        end
                        D(id) = D(id) + EL(n);
                    end
                end
            end
        end

        % Find total branch lengths of all descendant lineages originating
        % from a set of nodes (returns the TDD of node and all descendants,
        % stored in a vector D, which is indexed such that D(nodeIndex) is
        % the TDD of node(nodeIndex)

        function D = TotalDescendantDistance(obj,node,M)
            if isempty(obj.ancestors)
                error('Phylo object has no ancestral morphological data');
            end
            if nargin < 3
                [~,nV] = size(obj.ancestors);
                M = vertcat(obj.TaxonMatrix(nV),obj.ancestors);
            end
            [nR,~] = size(M);
            nT = length(obj.taxa);
            idx = node;
            D = zeros(nR,1);
            if idx > nT
                c = obj.nodes(idx).children;
                for i = 1:length(c)
                    D = D + obj.TotalDescendantDistance(c(i),M);
                    D(idx) = D(idx) + D(c(i)) + norm(M(idx,:)-M(c(i),:));
                end
            end
        end

        % find root node
        function R = FindRoot(obj)
            for i = 1:length(obj.nodes)
                if isempty(obj.nodes(i).parent)
                    R = i;
                end
            end
        end

        % interpolate a set of ancestor shapes along the ancestral lineage
        % of a node
        function S = InterpolateShapes(obj,node,BLs,M,EL)
            if isempty(obj.ancestors)
                error('Shape interpolation requires ancestral shape data');
            end
            if isempty(obj.edgeLengths)
                error('Shape interpolation requires edge length data');
            end

            N = length(BLs);
            if nargin < 4
                [~,nV] = size(obj.ancestors);
                M = vertcat(obj.TaxonMatrix(nV),obj.ancestors);
            end
            [~,nV] = size(M);
            S = zeros(N,nV);
            if nargin < 5
                EL = obj.NodeEL();
            end

            L = obj.Lineage(node);
            nL = length(L);
            T = zeros(nL,1);
            for i = 1:nL-1
                T(i+1) = T(i) + EL(L(i));
            end

            for i = 1:N
                if (BLs(i) - T(end)) > 0 && (BLs(i) - T(end)) < 0.0001
                    BLs(i) = T(end);
                end
                for j = 1:nL-1
                    id1 = L(j);
                    id2 = L(j+1);
                    if BLs(i) >= T(j) && BLs(i) <= T(j+1)
                        
                        d = BLs(i) - T(j);
                        r = d / (T(j+1)-T(j));
                        v = M(id2,:) - M(id1,:);
                        S(i,:) = M(id1,:) + r * v;
                    end
                end
            end
        end

        %% STATISTICS

        % pairwise convergence (default time informed (Grossnickle et al, 
        % preprint, doi: https://doi.org/10.1101/2022.10.18.512739))
        % can use legacy methodology from Stayton 2015
        function C = Convergence(obj,legacy)
            if isempty(obj.ancestors)
                error('Test of convergence requires ancestral morphological data');
            end
            if nargin < 2
                legacy = false;
            end
            if ~legacy && isempty(obj.edgeLengths)
                error('Time informed Convergence requires tree edge lengths. Use SetEdgeLengths or set legacy = true to use Stayton 2015 method');
            end

            % combined taxon ancestor matrix
            [~,nV] = size(obj.ancestors);
            M = vertcat(obj.TaxonMatrix(nV),obj.ancestors);
            % number of taxa
            nT = length(obj.taxa);
            % number of pairwise comparisons between taxa
            nP = nT * (nT-1) / 2;
            C = zeros(nP,4);

            % C3 denominator
            D3 = obj.EvoMorphDistance(1:nT);
            disp(D3);
            root = obj.FindRoot();
            % C4 denominator
            D4 = obj.TotalDescendantDistance(root,M);
            disp(D4);
            id = 0;

            L = cell(nT,1);
            for i = 1:nT
                L{i} = obj.Lineage(i);
            end

            % if not using legacy convergence, calculate total edge lengths
            if ~legacy
                EL = obj.NodeEL();
                T = cell(nT,1);
                for i = 1:nT
                    nL = length(L{i});
                    T{i} = zeros(nL,1);
                    for j = 1:nL-1
                        T{i}(j+1) = T{i}(j) + EL(L{i}(j));
                    end
                end
            end

            for i = 1:nT-1
                for j = i+1:nT
                    Dtip = norm(M(i,:)-M(j,:));
                    Dmax = 0;
                    [mrca,id1,id2] = obj.MRCA(L{i},L{j});
                    if legacy
                        for m = 1:id1
                            for n = 1:id2
                                if m ~= 1 && n ~= 1
                                    d = norm(M(L{i}(m),:) - M(L{j}(n),:));
                                    if d > Dmax
                                        Dmax = d;
                                    end
                                end
                            end
                        end
                    else
                        T1 = T{i}(id1) - (T{j}(id2) - T{j}(1:id2));
                        T2 = T{j}(id2) - (T{i}(id1) - T{i}(1:id1));
                        M1 = obj.InterpolateShapes(i,T1,M,EL);
                        M2 = obj.InterpolateShapes(j,T2,M,EL);

                        for m = 2:id1
                            d = norm(M(L{i}(m),:) - M2(m,:));
                            if T2(m) >= 0 && d > Dmax
                                Dmax = d;
                            end
                        end
                        for n = 2:id2
                            d = norm(M(L{j}(n),:) - M1(n,:));
                            if T1(n) >= 0 && d > Dmax
                                Dmax = d;
                            end
                        end
                    end

                    id = id+1;
                    if Dmax ==  0
                        C(id,:)  = [nan nan nan nan];
                    else
                        C(id,1) = 1 - (Dtip/Dmax);
                        C(id,2) = Dmax - Dtip;
                        C(id,3) = C(id,2) / D3(mrca);
                        C(id,4) = C(id,2) / D4(mrca);
                    end
                end
            end
        end

        %% PLOTS
        
        % plot tree
        function Plot(obj,minLength,showNames)
            if nargin < 3
                showNames = true;
            end
            [nR,~] = size(obj.edges);
            for i = 1:nR+1
                obj.nodes(i).Plot(obj.nodes,minLength);
                hold on
            end
            if showNames
                for i = 1:length(obj.taxa)
                    text(minLength * -0.5,i,obj.taxa(i).name, 'FontSize', 6);
                end
            end
            
            ax = gca;
            xL = ax.XLim;
            xL(2) = minLength * 3;
            ax.XLim = xL;
            
            axis off
        end
        
        % plot phylomorphospace
        function obj = Phylomorphospace(obj,TDS,ancestorFile,pcx,pcy)
            obj = obj.ImportEFAAncestors(ancestorFile);
            
            anc = obj.ancestors;
            
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
                if ~isempty(obj.nodes(i).parent)
                    plot([scr(i,pcx) scr(obj.nodes(i).parent,pcx)], [scr(i,pcy) scr(obj.nodes(i).parent,pcy)], 'k-');
                    hold on
                end
            end
            
            TDS.PlotClade(1,2);
            
        end
        

        % plot a tree with reconstructed ancestral shapes at each node
        function AncestorPhylo(obj,ancestorFile,minL,s)
            obj = obj.ImportEFAAncestors(ancestorFile);
            anc = obj.ancestors;
            
            [nr,~] = size(anc);
            
            for i = 1:nr
                ancShapes(i,1) = theoShapeN(1,anc(i,:));
            end
            
            for i = 1:length(obj.taxa)
                taxaShapes(i,1) = theoShapeN(1,obj.taxa(i).harm);
            end
            
            shapes = vertcat(taxaShapes,ancShapes);
            
            
            
            figure
            
            for i = 1:length(obj.nodes)
                obj.nodes(i).Plot(obj.nodes,minL);
                hold on
                x = obj.nodes(i).GetMaxBranchLength(obj.nodes) * minL * -1;
                [x,y] = shapes(i).draw(100,x,obj.nodes(i).y,minL*s);
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

        function PairwiseComparison(obj,Z,showNames)
            figure
            if nargin < 3
                showNames = true;
            end
            
            obj.Plot(2,showNames);
            hold on
            if showNames
                minX = 50;
            else
                minX = 1;
            end
            nT = length(obj.taxa);
            X = (1:nT) + minX; 
            Y = 1:nT;
            imagesc(X,Y,squareform(Z)');
            axis equal
            ax = gca;
            ax.XLim(2) = max(X);
        end

        function GroupComparison(obj,Z)
            groups = [obj.taxa.clade];
            nT = length(groups);
            g = cell(max(groups),1);
            for i = 1:nT
                if isempty(g{groups(i)})
                    g{groups(i)} = i;
                else
                    g{groups(i)} = horzcat(g{groups(i)},i);
                end
            end

            X = 1:nT; 
            Y = 1:nT;

            nG = length(g);
            xgrids = cell(nG,nG);
            ygrids = cell(nG,nG);
            zgrids = cell(nG,nG);

            ids = zeros(nG+1,1);
            
            for i = 1:nG
                ids(i+1) = ids(i) + length(g{i});
            end
            gap = 2;
            for i = 1:nG
                for j = 1:nG
                    xgrids{i,j} = X(ids(i)+1:ids(i+1)) + gap * (i-1);
                    ygrids{i,j} = Y(ids(j)+1:ids(j+1)) + gap * (j-1);
                    zgrids{i,j} = zeros(length(g{i}),length(g{j}));
                    for m = 1:length(g{i})
                        mID = g{i}(m);
                        for n = 1:length(g{j})
                            nID = g{j}(n);
                            if mID == nID
                                zgrids{i,j}(m,n) = 0;
                            else
                                if mID <= nID
                                    I = mID;
                                    J = nID;
                                else
                                    I = nID;
                                    J = mID;
                                end
    
                                ID = (I - 1) * (nT - I/2) + J - I;
                                zgrids{i,j}(m,n) = Z(ID);
                            end
                        end
                    end
                end
            end

            figure

            for i = 1:nG
                for j = 1:nG
                    imagesc(xgrids{i,j},ygrids{i,j},zgrids{i,j}');
                    hold on
                end
            end
            axis equal
            xlim([min(xgrids{1,1}) max(xgrids{nG,nG})]);
            ylim([min(ygrids{1,1}) max(ygrids{nG,nG})]);

            axis off
            
        end
    end
end

