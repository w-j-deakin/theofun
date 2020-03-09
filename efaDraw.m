function [x,y] = efaDraw(harm, n)
% Draw an outline with n co-ordinates based on harm efa harmonic data

[nHarm,~] = size(harm);

for i = 1:n
    for j = 1:nHarm
        X(i,j) = (harm(j,1) * cos(j * (i*2*pi)/n)) + (harm(j,2) * sin(j * (i*2*pi)/n));
        Y(i,j) = (harm(j,3) * cos(j * (i*2*pi)/n)) + (harm(j,4) * sin(j * (i*2*pi)/n));
    end 
    x(i) = sum(X(i,:));
    y(i) = sum(Y(i,:));
end

end

