function [isIntersect] = lineIntersect(line1,line2)


% Line segments intersect if:
%       - Bounding boxes intersect
%       - line a intersects with segment b
%       - line b intersects with segment a


aX1 = min(line1(:,1));
aX2 = max(line1(:,1));

aY1 = min(line1(:,2));
aY2 = max(line1(:,2));

bX1 = min(line2(:,1));
bX2 = max(line2(:,1));

bY1 = min(line2(:,2));
bY2 = max(line2(:,2));



if aX1 <= bX2 && bX1 <= aX2 && aY1 <= bY2 && bY1 <= aY2
    
    aS  = [(line1(1,1) - line1(2,1)) (line1(1,2) - line1(2,2))];
    aSb = [(line2(1,1) - line1(2,1)) (line2(1,2) - line1(2,2));
           (line2(2,1) - line1(2,1)) (line2(2,2) - line1(2,2))];
    bS  = [(line2(1,1) - line2(2,1)) (line2(1,2) - line2(2,2))];
    bSa = [(line1(1,1) - line2(2,1)) (line1(1,2) - line2(2,2));
           (line1(2,1) - line2(2,1)) (line1(2,2) - line2(2,2))];
    
    cpA1 = (aS(1) * aSb(1,2)) - (aS(2) * aSb(1,1));
    cpA2 = (aS(1) * aSb(2,2)) - (aS(2) * aSb(2,1));
    cpB1 = (bS(1) * bSa(1,2)) - (bS(2) * bSa(1,1));
    cpB2 = (bS(1) * bSa(2,2)) - (bS(2) * bSa(2,1));
    
    if (cpA1 > 0 && cpA2 < 0)||(cpA1 < 0 && cpA2 > 0)
        if (cpB1 > 0 && cpB2 < 0)||(cpB1 < 0 && cpB2 > 0)
            isIntersect = true;
        else
            isIntersect = false;
        end
    else
        isIntersect = false;
    end
    
else
    isIntersect = false;
end

end

