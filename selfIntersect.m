function [isIntersect] = selfIntersect(harm)

[x,y] = efaDraw(harm,200);

inComp = true;
i = 1;

while inComp
    if i == 200
        isIntersect = false;
        inComp = false;
    else
        lineA = [x(i) y(i);
                 x(i+1) y(i+1)];
        j = i + 1;
        run = true;
        while run
            if j == 200
                lineB = [x(j) y(j);
                         x(1) y(1)];
                if lineIntersect(lineA,lineB)
                    isIntersect = true;
                    inComp = false;
                    run = false;
                else
                    run = false;
                end
            else
                lineB = [x(j) y(j);
                         x(j+1) y(j+1)];
                if lineIntersect(lineA,lineB)
                    isIntersect = true;
                    inComp = false;
                    run = false;
                else
                    j = j + 1;
                end
            end    
        end
        i = i + 1;
    end
end


end

