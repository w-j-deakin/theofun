classdef mesh2D
    
    properties
        x
        y
        outline
        ele = ele2D.empty(0,1)
        nEle
    end
    
    methods
        %% CONSTRUCTOR
        function obj = mesh2D(inX,inY,inEle,OL)
            
            if nargin == 0
                obj.x = nan;
                obj.y = nan;
                obj.outline = nan;
                obj.nEle = nan;
            else
                [rX,cX] = size(inX);
                if rX > 1 && cX > 1
                    error('input x must be a vector');
                elseif cX > 1
                    inX = inX';
                end  

                [rY,cY] = size(inY);
                if rY > 1 && cY > 1
                    error('input y must be a vector');
                elseif cY > 1
                    inY = inY';
                end 
                if length(inY) ~= length(inX)
                    error('Vectors x and y must be the same length');
                end

                if inEle(1).test_ele2D
                    obj.ele = inEle;
                end

                [rOL,cOL] = size(OL);
                if rOL > 1 && cOL > 1
                    error('Input outline must be an integer vector');
                elseif cOL > 1
                    OL = OL';
                end

                if length(OL) > length(inX)
                    error('Outline vector must be shorter than x and y vectors');
                elseif ~isinteger(OL)
                    error('Outline vector must have integer elements');
                end

                oX = inX(OL);
                oY = inY(OL);

                area = polyarea(oX,oY);
                inX = inX/sqrt(area);
                inY = inY/sqrt(area);
                oX = inX(OL);
                oY = inY(OL);
                poly = polyshape(oX,oY);
                [centX, centY] = centroid(poly);
                inX = inX - centX;
                inY = inY - centY;

                obj.x = inX;
                obj.y = inY;
                obj.outline = OL;
                nEle = length(inEle);
                obj.nEle = nEle;
            end
        end
        
        %% PLOTS
        % Draw mesh
        
        function Draw(obj,col,pos,scale)
            if nargin < 2
                col = ones(obj.nEle,3);
                col = col.*0.75;
            end
            if nargin < 3
                pos = [0 0];
            end
            if nargin < 4
                scale = 1;
            end
           
            for i = 1:obj.nEle
                fill(obj.x(obj.ele(i).nodes)*scale+pos(1),obj.y(obj.ele(i).nodes)*scale+pos(2),col(i,:),'LineStyle','none');
                hold on
            end
%             fill(obj.x(obj.outline),obj.y(obj.outline),[0.75 0.75 0.75],'LineStyle','none'),
        end
        
        % Draw FE Model
        
        function DrawFEM(obj,pos,scale)
            
            if nargin < 2
                pos = [0 0];
            end
            if nargin < 3
                scale = 1;
            end
            vms = obj.FEA;
            disp('vms calculated');
            colormap jet
            cmap = colormap;
            [nC,~] = size(cmap);
            
            mV = prctile(vms,99);
            
            R = mV - min(vms);
            
            r = nC/R;
            
            col = zeros(obj.nEle,3);
            
            for i = 1:obj.nEle
                cI = floor((vms(i) - min(vms)) * r) + 1;
                if cI > nC
                    col(i,:) = [0.75 0.75 0.75];
                else
                    col(i,:) = cmap(cI,:);
                end
            end
            
            obj.Draw(col,pos,scale);
        end
        
        %% AUTOMATIC MODELS
        % Auto constraints
        
        function [c,nrml] = autoCon(obj,range)
            
            if nargin < 2
                range = [ 0.25 -0.85;
                         -0.25  0.85 ];
            end
            
            nOL = length(obj.outline);
            
            oX = obj.x(obj.outline);
            oY = obj.y(obj.outline);
            clockwise = ispolycw(oX,oY);

            nrml = zeros(nOL,2);
            nNrml = zeros(nOL,1);
            minX = 0;
            maxX = 0;

            % Automatic constraints
            for i = 1:nOL
                % Define vectors
                if i == 1
                    A(1) = oX(nOL) - oX(1);
                    A(2) = oY(nOL) - oY(1);
                    B(1) = oX(2) - oX(1);
                    B(2) = oY(2) - oY(1);
                elseif i == nOL
                    A(1) = oX(nOL-1) - oX(nOL);
                    A(2) = oY(nOL-1) - oY(nOL);
                    B(1) = oX(1) - oX(nOL);
                    B(2) = oY(1) - oY(nOL);
                else
                    A(1) = oX(i-1) - oX(i);
                    A(2) = oY(i-1) - oY(i);
                    B(1) = oX(i+1) - oX(i);
                    B(2) = oY(i+1) - oY(i);
                end

                % Bisector of two vectors A and B = |A|*B + |B|*A
                nrml(i,:) = (norm(A)*B + norm(B)*A);
                nNrml(i) = norm(nrml(i,:));


                aveX = (A(1) - B(1))/2;
                aveY = (A(2) - B(2))/2;

                if(~clockwise)
                    aveX = -aveX;
                    aveY = -aveY;
                end

                if aveX <= 0 && nrml(i,2) <= 0
                    nrml(i,:) = nrml(i,:) * -1;
                end

                if aveY <= 0 && nrml(i,1) >= 0
                    nrml(i,:) = nrml(i,:) * -1;
                end

                if aveX >= 0 && nrml(i,2) >= 0
                    nrml(i,:) = nrml(i,:) * -1;
                end

                if aveY >= 0 && nrml(i,1) <= 0
                    nrml(i,:) = nrml(i,:) * -1;
                end

                nrml(i,:) = nrml(i,:) / nNrml(i);

                if oX(i) < minX && nrml(i,2) > 0 && nrml(i,1) > range(1,2) && nrml(i,1) < range(1,1)
                    minX = oX(i);
                    c(1) = obj.outline(i);     
                elseif oX(i) > maxX && nrml(i,2) > 0 && nrml(i,1) < range(2,2) && nrml(i,1) > range(2,1) 
                    maxX = oX(i);
                    c(2) = obj.outline(i);
                end
            end
        end
        
        % Auto force
        
        function [f,c] = autoF(obj,c,nrml)
            if nargin < 2
                [c,nrml] = obj.autoCon;
            elseif nargin < 3
                [~,nrml] = obj.autoCon;
            end
            nOL = length(obj.outline);
            
            nX = obj.x;
            nY = obj.y;
            
            oX = nX(obj.outline);
            oY = nY(obj.outline);
            
            
            lngth(1) = (nX(c(2)) - nX(c(1))) * 0.25;
            lngth(2) = (nY(c(2)) - nY(c(1))) * 0.25;

            fP(1) = nX(c(1)) + lngth(1);
            fP(2) = nX(c(1)) + lngth(2);

            fLU(1) = - (lngth(2)/lngth(1))*1000;
            fLU(2) = 1000;

            fLD(1) = (lngth(2)/lngth(1)) * 1000;
            fLD(2) = - 1000;

            [xUi,yUi] = polyxpoly([fP(1) fLU(1)],[fP(2) fLU(2)],oX,oY);
            [xDi,yDi] = polyxpoly([fP(1) fLD(1)],[fP(2) fLD(2)],oX,oY);

            [nU,~] = size(xUi);
            [nD,~] = size(xDi);

            if nU ~= 0
                iU = zeros(nU,1);
            end
            if nD ~= 0
                iD = zeros(nD,1);
            end


            fN = 0;

            for i = 1:nD
                minDist = 1;
                for j = 1:nOL
                    dist = norm([xDi(i) yDi(i)] - [oX(j) oY(j)]);
                    if dist < minDist
                        minDist = dist;
                        iD(i) = j;
                    end
                end
                if nrml(iD(i),2) > 0
                    fN = iD(i);
                end
            end
            for i = 1:nU
                minDist = 1;
                for j = 1:nOL
                    dist = norm([xUi(i) yUi(i)] - [oX(j) oY(j)]);
                    if dist < minDist
                        minDist = dist;
                        iU(i) = j;
                    end
                end
                if nrml(iU(i),2) > 0
                    fN = iU(i);
                end
            end

            if fN == 0
                f = [0 0 0];
            else
                fI = double(obj.outline(fN));
                fVec = nrml(fN,:)/norm(nrml(fN,:));
                f = [fI fVec(1) fVec(2)];
            end
        end
        
        % Create Node-Constraint Vector
        
        function C = vecCon(obj,c)
            nN = length(obj.x);
            C = zeros(nN,1);
            C(c) = 1;
        end
        
        % Create Node-Force Vector
        
        function F = vecF(obj,f)
            if f(1) == 0
                F = nan;
            else
                nN = length(obj.x);
                F = zeros(nN,2);
                F(f(1),:) = [f(2) f(3)];
            end
        end
        
        % Random sample of auto constraints
        
        function cV = randCon(obj,n,r)
            if nargin < 3
                r = 0.1;
            end
            
            N = length(obj.outline);
            rn = r*N;
            if mod(round(rn),2)
                nP = (round(rn)-1)/2;
            else
                nP = round(rn)/2;
            end
            
            midC = obj.autoCon([0.25 -0.75;
                               -0.25  0.75]);
            jntC = (midC(1) - nP):(midC(1) + nP);
            tthC = (midC(2) - nP):(midC(2) + nP);
            
            for i = 1:2*nP+1
                if jntC(i) < 1
                    jntC(i) = N + jntC(i);
                elseif jntC(i) > N
                    jntC(i) = jntC(i) - N;
                end
                
                if tthC(i) < 1
                    tthC(i) = N + tthC(i);
                elseif tthC(i) > N
                    tthC(i) = tthC(i) - N;
                end 
            end
            
            cV = zeros(n,2);
            for i = 1:n
                cV(i,1) = jntC(randi(2*nP+1));
                cV(i,2) = tthC(randi(2*nP+1));
            end
        end
        
        % Random sample auto force
        
        function fV = randForce(obj,n,r,aR)
            if nargin < 3
                r = 0.1;
            end
            if nargin < 4
                aR = [-pi/4 pi/4];
            end
            
            ar = aR(2) - aR(1);
            
            f = obj.autoF;
            midN = f(1);
            
            N = length(obj.outline);
            rn = r*N;
            if mod(round(rn),2)
                nP = (round(rn)-1)/2;
            else
                nP = round(rn)/2;
            end
            
            nV = (midN - nP):(midN + nP);
            
            for i = 1:2*nP+1
                if nV(i) < 1
                    nV(i) = N + nV(i);
                    disp(nV(i));
                elseif nV(i) > N
                    nV(i) = nV(i) - N;
                    disp(nV(i));
                end
            end
            
            fV = zeros(n,3);
            for i = 1:n
                fV(i,1) = nV(randi(2*nP+1));
                A = rand*ar + aR(1);
                V = Rotate(f(2:3),A);
%                 fV(i,2) = sin(A);
%                 fV(i,3) = cos(A);
                fV(i,2) = V(1);
                fV(i,3) = V(2);
            end
        end
        
        % Viusalise default constraints
        
        function DrawDefaultInputs(obj,randCon,randForce,ax)
            if nargin < 2
                randCon = true;
            end
            if nargin < 3
                randForce  = true;
            end
            if nargin < 4
                ax = axes;
            end
            axes(ax);
            obj.Draw;
            hold on
            plot(vertcat(obj.x(obj.outline), obj.x(obj.outline(1))),vertcat(obj.y(obj.outline), obj.y(obj.outline(1))),'k-','LineWidth',0.5);
            
            [~,nrml] = obj.autoCon;
            if randCon
                cV = obj.randCon(1000);
            else
                cV = obj.autoCon;
            end
            fV = obj.randForce(1000);
            
            J = unique(cV(:,1));
            T = unique(cV(:,2));
            F = unique(fV(:,1));
            fP = sort(F);
            
            
            if randForce
                for i = 1:length(F)
                    X = obj.x(F(i));
                    Y = obj.y(F(i));

%                     X1 = nrml(F(i),1) + sin(-pi/4);
%                     X2 = nrml(F(i),1) + sin(pi/4);
                     
%                     Y1 = nrml(F(i),2) + cos(-pi/4);
%                     Y2 = nrml(F(i),2) + cos(pi/4);
                    v1 = Rotate(nrml(F(i),:),-pi/4);
                    v2 = Rotate(nrml(F(i),:), pi/4);

                    V1 = 0.5*v1 / norm(v1);
                    V2 = 0.5*v2 / norm(v2);
                    
                    fX(i,:) = [X (X+V1(1)) (X+0.5*nrml(F(i),1)) (X+V2(1))];
                    fY(i,:) = [Y (Y+V1(2)) (Y+0.5*nrml(F(i),2)) (Y+V2(2))];

                    fill(fX(i,:),fY(i,:),hsv2rgb([0 0.35 1]),'LineStyle','none');
                end
%                 for i = 1:length(F)
% %                     fill(fX(i,:),fY(i,:),'w','FaceAlpha',0,'EdgeColor',hsv2rgb([0 0.75 1]));
%                 end
                plot(obj.x(fP),obj.y(fP),'r-','LineWidth',2);
                
                for i = 1:length(J)
                    X = obj.x(J(i));
                    Y = obj.y(J(i));
                    plot(X,Y,'ko','MarkerFaceColor','b');
                end

                for i = 1:length(T)
                    X = obj.x(T(i));
                    Y = obj.y(T(i));
                    plot(X,Y,'ko','MarkerFaceColor','b');
                end
            end
            
            jP = 0;
            tP = 0;
            for i = 1:length(obj.outline)
                for j = 1:length(J)
                    if J(j) == i
                        jP = horzcat(jP,i);
                    end
                end
                
                for j = 1:length(T)
                    if T(j) == i
                        tP = horzcat(tP,i);
                    end
                end
            end
            jP(1) = [];
            tP(1) = [];
            
            for i = 2:length(jP)
                if jP(i) - jP(i-1) ~= 1
                    jP = horzcat(jP(i:end),jP(1:i-1));
                end
            end
            for i = 2:length(tP)
                if tP(i) - tP(i-1) ~= 1
                    tP = horzcat(tP(i:end),tP(1:i-1));
                end
            end
            plot(obj.x(tP),obj.y(tP),'b-','LineWidth',2);
            plot(obj.x(jP),obj.y(jP),'b-','LineWidth',2);
            axis off
        end
        
        %% FEA
        % Global stiffness matrix
        
        function [dK,D,B,shft,K,C] = gSM(obj,t,E,v,C)
            
            if nargin < 2
                t = 0.01;
            end
            if nargin < 3
                E = 2*(10^9);
            end
            if nargin < 4
                v = 0.3;
            end
            if nargin < 5
                c = obj.autoCon;
                C = obj.vecCon(c);
            end
            
            nN = length(obj.x);
            nE = obj.nEle;
            

            nC = 0;
            shft = zeros(nN,1);

            for i = 1:nN
                shft(i) = i - nC;
                if(C(i)) == 1
                    nC = nC + 1;
                end
            end

            % D is a matrix used in the stiffness matrix calculation, which inputs the
            % youngs modusu and poisson ratio. It characterises the relationship
            % between stress and strain in 2D, such that:
            %
            %     stress vector = D * strain vector

            a = E / (1 - v*v);
            b = (1 - v) / 2;

            D = a * [ 1 v 0 ;
                      v 1 0 ;
                      0 0 b ];

            % K will be the global stiffness matrix
            % For speed, we will preallocate it with constraints emitted
            nK = (nN - nC) * 2;
            K = zeros(nK,nK);
            ke = zeros(6,6,nE);
            B = zeros(3,6,nE);

            % Loop through every element
            for i = 1:nE
                % Get element points
                idx(1) = obj.ele(i).nodes(1);
                idx(2) = obj.ele(i).nodes(2);
                idx(3) = obj.ele(i).nodes(3);

                x1 = obj.x(idx(1));
                y1 = obj.y(idx(1));

                x2 = obj.x(idx(2));
                y2 = obj.y(idx(2));

                x3 = obj.x(idx(3));
                y3 = obj.y(idx(3));


                % The stiffness matrix for each element (ke) can be calculated as
                % follows:
                %
                % ke = Te * Ae * Bt * D * B
                %
                % Where: 
                % Te = element thickness (function input argument)
                % Ae = element area
                % D = material matrix already calculated
                % B = the relationship between nodal displacement and strain
                % Bt = B transposed

                % Calculate Area = 0.5 * |detJ| (determinant of matrix jacobian)
                % For simplicity, let x12 = x1 - x2, y21 = y2 - y1, etc.
                % in 2D triangular examples, detJ = x13*y23 - y13*x23
                x13 = x1 - x3;
                x21 = x2 - x1;
                x23 = x2 - x3;
                x32 = x3 - x2;

                y12 = y1 - y2;
                y13 = y1 - y3;
                y23 = y2 - y3;
                y31 = y3 - y1;

                detJ = (x13*y23) - (y13*x23);

                Ae = 0.5 * abs(detJ);

                % B is defined as:

                B(:,:,i) = (1/detJ) * [ y23   0   y31   0   y12   0  ; 
                                  0   x32   0   x13   0   x21 ;
                                 x32  y23  x13  y31  x21  y12 ];

                temp = B(:,:,i);               
                Bt = transpose(temp);

                % Perform matrix calculations in correct order first
                temp1 = Bt * D;
                temp2 = temp1 * temp;
                ke(:,:,i) = t * Ae * temp2;

                % Find shifted indices
                for j = 1:3
                    shift(2*j - 1) = 2*shft(idx(j)) - 1;
                    shift(2*j) = 2*shft(idx(j));
                end

                % Don't add to global stiffness matrix if the node is a constraint
                for j = 1:6
                    for J = 1:6
                        if C(idx(ceil(j/2))) == 0 && C(idx(ceil(J/2))) == 0
                            K(shift(j),shift(J)) = K(shift(j),shift(J)) + ke(j,J,i);
                        end
                    end
                end
            end
            
            dK = decomposition(K);
            
        end
        
        % FEA from stiffness matrix
        
        function vms = feaK(obj,F,C,dK,D,B,shft)
            nN = length(obj.x);
            nE = obj.nEle;
            nK = 2*(nN-2);

            % FORCES

            % Define force vector
            f = zeros(nK,1);
            for i = 1:nN
                if C(i) == 0
                    I = shft(i);
                    f(2*I - 1) = F(i,1);
                    f(2*I) = F(i,2);
                end
            end

            % LINEAR EQUATION

            % Solve eaquation for displacements (Q)

            Q = dK\f;

            if Q(1) == 0 || isnan(Q(1))
                vms = 0;
                return
            end

            % DISPLACEMENTS

            % Add displacements to nodes
            add = zeros(2,1);
            for i = 1:nN
                if C(i) == 1
                    Q = cat(1,Q(1:(2*(i-1))),add,Q(((2*i)-1):end));
                end
            end

            % STRESS AND STRAIN

            % Back calculate stress and strain from displacements
            sig = zeros(nE,3);
            vms = zeros(nE,1);

            for i = 1:nE
                q1 = (obj.ele(i).nodes(1))*2 - 1;
                q2 = (obj.ele(i).nodes(1))*2;
                q3 = (obj.ele(i).nodes(2))*2 - 1;
                q4 = (obj.ele(i).nodes(2))*2;
                q5 = (obj.ele(i).nodes(3))*2 - 1;
                q6 = (obj.ele(i).nodes(3))*2;

                % Use B and D to calculate strain, tehn stress vectors
                temp1 = B(:,:,i);
                temp2 = D*temp1;
                temp3 = [Q(q1);
                         Q(q2);
                         Q(q3);
                         Q(q4);
                         Q(q5);
                         Q(q6)];

                % Calculate stress vector
                sig(i,:) = temp2*temp3;

                %VMS formula
                vms(i) = ((sig(i,1)^2) + (sig(i,2)^2) - (sig(i,1)*sig(i,2)) + (3*(sig(i,3)^2)))^0.5;
            end
        end
        
        % Automatic FEA, start to finish.
        
        function vms = FEA(obj,t,E,v,c,f)
            
            if nargin < 2
                t = 0.01;
            end
            if nargin < 3
                E = 2*(10^9);
            end
            if nargin < 4
                v = 0.3;
            end
            if nargin < 6
                if nargin < 5
                    [ct,nrml] = obj.autoCon;
                    if length(ct) == 2 && ct(1) > 0 && ct(2) > 0
                        c = ct;
                    else
                        vms = zeros(obj.nEle,1);
                        return
                    end
                    f = obj.autoF(c,nrml);
                else
                    f = obj.autoF(c);
                end
            end
            
            C = obj.vecCon(c);
            F = obj.vecF(f);
            
            if isnan((F))
                vms = zeros(obj.nEle,1);
            else
                [dK,D,B,shft] = obj.gSM(t,E,v,C);
                vms = obj.feaK(F,C,dK,D,B,shft);
            end
        end
        
        % Bootstrap FEA
        
        function vms = btstrpFEA(obj,n,N,r,aR,t,E,v)
            
            if nargin < 3
                N = [1 1];
            end
            if nargin < 4
                r = 0.1;
            end
            if nargin < 5
                aR = [-pi/4 pi/4];
            end
            if nargin < 6
                t = 0.01;
            end
            if nargin < 7
                E = 2*(10^9);
            end
            if nargin < 8
                v = 0.3;
            end
            
            [dK,D,B,shft,~,C] = obj.gSM(t,E,v);
            fV = obj.randForce(n,r,aR);
            vms = zeros(obj.nEle,n);
            p = zeros(n,1);
            S = ['Mesh ',num2str(N(1)),'/',num2str(N(2))];
            disp(S);
            disp('VMS Bootstrap completion:');
            for i = 1:n
                F = obj.vecF(fV(i,:));
                vms(:,i) = obj.feaK(F,C,dK,D,B,shft); 
                p(i) = i*100/n;
                if (p(i) >= 25 && p(i-1) < 25) || (p(i) >= 50 && p(i-1) < 50) || (p(i) >= 75 && p(i-1) < 75)
                    S = [num2str(p(i)),'%'];
                    disp(S);
                end
            end
            clc;
        end
        
        %% ROTATIONAL INERTIA
        
        % Find element mass and center
        function [M,C] = eleMC(obj)
            nE = obj.nEle;
            
            C = zeros(nE,2);
            M = zeros(nE,1);
            
            for i = 1:nE
                % Get element points
                idx(1) = obj.ele(i).nodes(1);
                idx(2) = obj.ele(i).nodes(2);
                idx(3) = obj.ele(i).nodes(3);

                x1 = obj.x(idx(1));
                y1 = obj.y(idx(1));

                x2 = obj.x(idx(2));
                y2 = obj.y(idx(2));

                x3 = obj.x(idx(3));
                y3 = obj.y(idx(3));

                % Find element center
                C(i,1) = (x1 + x2 + x3) / 3;
                C(i,2) = (y1 + y2 + y3) / 3;

                % Find area of element
                M(i) = abs(x1*(y2-y3) + x2*(y3-y1) + x3*(y1-y2));
            end
        end
        
        % RI, v using M & C
        function [v,RI] = MCRI(obj,M,C,p0,p)
            nE = obj.nEle;
            dist = zeros(nE,1);
            rie = zeros(nE,1);
            for i = 1:nE
                dist(i) = norm(p0 - C(i,:));
                rie(i) = M(i) * dist(i)^2;
            end
            RI = sum(rie);
            L = norm(p - p0);
            v = sqrt(2)*L/sqrt(RI);
        end
        
        % Calculate RI
        
        function RI = RI(obj,p0)
            
            if nargin < 2
                c = obj.autoCon;
                p0 = [obj.x(c(1)) obj.y(c(1))];
            end
            
            nE = obj.nEle;
            
            dist = zeros(nE,1);
            area = zeros(nE,1);
            riE = zeros(nE,1);
            
            for i = 1:nE

                % Get element points
                idx(1) = obj.ele(i).nodes(1);
                idx(2) = obj.ele(i).nodes(2);
                idx(3) = obj.ele(i).nodes(3);

                x1 = obj.x(idx(1));
                y1 = obj.y(idx(1));

                x2 = obj.x(idx(2));
                y2 = obj.y(idx(2));

                x3 = obj.x(idx(3));
                y3 = obj.y(idx(3));

                % Find element center
                cnt(1) = (x1 + x2 + x3) / 3;
                cnt(2) = (y1 + y2 + y3) / 3;

                % Find distance from joint
                dist(i) = norm(cnt - p0);

                % Find area of element
                area(i) = abs(x1*(y2-y3) + x2*(y3-y1) + x3*(y1-y2));

                % Find RI of element
                riE(i) = area(i) * dist(i)^2;
            end
            % RI = sum of RI for all elements
            RI = sum(riE);
        end
        
        % Calculate speed per unit energy (rotational efficiency)
        
        function [v,RI] = rotV(obj,p0,p)
            if nargin < 3
                ct = obj.autoCon;
                if length(ct) == 2 && ct(1) > 0 && ct(2) > 0
                    c = ct;
                else
                    v = 0;
                    RI = 0;
                    return
                end
                p = [obj.x(c(2)) obj.y(c(2))];
                if nargin < 2
                    p0 = [obj.x(c(1)) obj.y(c(1))];
                end
            end
            RI = obj.RI(p0);
            L = norm(p - p0);
            v = sqrt(2)*L/sqrt(RI);
            
        end
        
        % Bootstrap RI
        
        function [v,RI] = btstrpRI(obj,n,N,r)
            if nargin < 2
                N = [1 1];
            end
            if nargin < 3
                r = 0.1;
            end
            
            cV = obj.randCon(n,r);
            v = zeros(n,1);
            RI = zeros(n,1);
            [M,C] = obj.eleMC;
            p = zeros(n,1);
            S = ['Mesh ',num2str(N(1)),'/',num2str(N(2))];
            disp(S);
            disp('Rotational Inertia Bootstrap completion:');
            for i = 1:n
                p0 = [obj.x(cV(i,1)) obj.y(cV(i,1))];
                p1 = [obj.x(cV(i,2)) obj.y(cV(i,2))];
                
                [v(i),RI(i)] = obj.MCRI(M,C,p0,p1);
                
                p(i) = i*100/n;
                if (p(i) >= 25 && p(i-1) < 25) || (p(i) >= 50 && p(i-1) < 50) || (p(i) >= 75 && p(i-1) < 75)
                    S = [num2str(p(i)),'%'];
                    disp(S);
                end
            end
            clc;
        end
        
        %%
    end
end

