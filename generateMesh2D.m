function mesh = generateMesh2D(harm,nEle)
% Draw test outline
[xS , yS] = efaDraw(harm,100);
sample = polyshape(xS,yS);
[centX, centY] = centroid(sample);

% Calculate number of points based on number of elements
area = polyarea(xS,yS);
perim = perimeter(sample);

A = area/nEle;
dsq = (2*A)/(sin(pi/3));
d = sqrt(dsq);

n = round(perim/d);

% Update d
d = perim/n;
h = d*sin(pi/3);

% Draw model outline
[x , y] = efaDraw(harm,n);

for i = 1:n
    x(i) = x(i) - centX;
    y(i) = y(i) - centY;
end

points = horzcat(x',y');


%% Define triangle point grid
k = 0;
K = 0;
No = 1;
row = 1;
xmax = 1.1*max(x);
xmin = 1.1*min(x);
ymax = 1.1*max(y);
ymin = 1.1*min(y);

while k == 0
    if mod(row,2) == 1
         if (xmin+(d/2))+(K*d) < xmax
            grid(No,:) = [(xmin+(d/2))+(K*d) ymin];
            No = No+1;
            K = K+1;
         else
            ymin = ymin + h;
            K = 0;
            row = row + 1;
            if ymin > ymax
                k = 1;
            end
        end
    else
        if xmin+(K*d) < xmax
            grid(No,:) = [xmin+(K*d) ymin];
            No = No+1;
            K = K+1;
         else
            ymin = ymin + h;
            K = 0;
            row = row + 1;
            if ymin > ymax
                k = 1;
            end
        end
    end
end
%%
%% Define Mesh

% Find all points in grid that lie inside shape
in = inpolygon(grid(:,1),grid(:,2),x,y);

% If any points are closer than d/1.75 to a point on the outline, ignore them
k = 0;
for i = 1:length(in)
    if in(i) == 1
        L = linspace(0,2*pi,10);
        xv = (d/(1.75))*(cos(L)')+grid(i,1);
        yv = (d/(1.75))*(sin(L)')+grid(i,2);
        close = inpolygon(x,y,xv,yv);
        if any(close) == 0
            k = k+1;
            inNodes(k,:) = [grid(i,1) grid(i,2)];
        end
    end     
end

% Node matrix = outline points + inside points
nodes = vertcat(points,inNodes);

% Delaunay triangulation to generate elements
ele = delaunay(nodes);
[eleNo,~] = size(ele);

% Erroneous elements outside original shape deleted, and generate element
% object array
mid = zeros(eleNo,2);
Ele2D = ele2D.empty(0,1);

k = 0;
for i = 1:eleNo
    x1 = nodes(ele(i,1),1);
    x2 = nodes(ele(i,2),1);
    x3 = nodes(ele(i,3),1);
    y1 = nodes(ele(i,1),2);
    y2 = nodes(ele(i,2),2);
    y3 = nodes(ele(i,3),2);
    mid(i,:) = [((x1+x2+x3)/3)  ((y1+y2+y3)/3)];
    if inpolygon(mid(i,1),mid(i,2),x,y) == 1
        k = k+1;
        Ele2D(k) = ele2D(ele(i,:));
    end
end

mesh = mesh2D(nodes(:,1),nodes(:,2),Ele2D,int16(1:n));


end

