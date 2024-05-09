function [N] = EFANormal(harm,t)

nt = length(t);
[nh,~] = size(harm);

N = zeros(nt,2);

for i = 1:nt
    xpt = 0;
    ypt = 0;
    for j = 1:nh
        xpt = xpt + j * (harm(j,2) * cos(j*t(i)) - harm(j,1) * sin(j*t(i)));
        ypt = ypt + j * (harm(j,4) * cos(j*t(i)) - harm(j,3) * sin(j*t(i)));
    end

    alpha = atan2(ypt,xpt);
    N(i,:) = [cos(alpha-pi/2) sin(alpha-pi/2)];
end


end