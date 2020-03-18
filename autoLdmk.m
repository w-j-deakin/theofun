function [ldmk] = autoLdmk(image,n)


% Load image
img = imread(image);

[imgH,imgW] = size(img);

% Starting from the top left corner, move diagonally until a pixel is found
xS = 2;
yS = 2;
k = 0;
isSearching = true;

while(isSearching)
    if xS >= imgH
        k = k + 1;
        xS = 2 + k;
    end
    if yS >= imgW
        k = k + 1;
        yS = 2 + k;
    end
    if img(xS,yS) == 0
        isSearching = false;
    else
        xS = xS + 1;
        yS = yS + 1;
    end
end

% Use this pixel as a starting point to ldmk the outline
isComplete = false;
pix(1,:) = [xS yS];
idx = 1;
vis = zeros(imgW,imgH);
mina = pi/4;

L = linspace(0,2*pi,10);
xv = 2*(cos(L)')+xS;
yv = 2*(sin(L)')+yS;


while(~isComplete)
    % Set last pixel as visited
    vis(pix(idx,1),pix(idx,2)) = 1;
    
    isFound = false;
    
    while ~isFound    
        % Use spiral to find next pixel
        pix(idx+1,:) = spiral(pix(idx,:),img,vis);
        
        % Check if neighbour pixel is not behind last pixel (after 5
        % pixels)
        if idx > 5
            p5 = pix(idx-5,:);
                
            P0 = pix(idx,:) - p5;
            P1 = pix(idx,:) - pix(idx+1,:);
            [~,alpha] = findAngle2D(P0,P1);
            
            if alpha < mina
                vis(pix(idx+1,1),pix(idx+1,2)) = 1;
            else
                isFound = true;
            end
            
        else
            isFound = true;
        end
    end 
    idx = idx + 1;
    
    % Check if complete
    close = inpolygon(pix(idx,1),pix(idx,2),xv,yv);
    if idx > 5 && close == 1
        isComplete = true;
    end
end

% Downsample co-prdinates to n
[nPix,~] = size(pix);

res = nPix/n;

for i = 1:n
    idxL = round((i-1)*res) + 1;
    ldmk(i,:) = [pix(idxL,1) -pix(idxL,2)];
end

if ispolycw(ldmk(:,1),ldmk(:,2))
    ldmk = flipud(ldmk);
end

end

