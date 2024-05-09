function [at] = EFASecondDerivative(harm, t)

[nH,~] = size(harm);
nT = length(t);

at = zeros(nT,1);

for i = 1:nT
    ft = 0;
    gt = 0;
    fpt = 0;
    gpt = 0;
    for j = 1:nH
        ft = ft + j * (harm(j,4)*cos(j*t(i)) - harm(j,3)*sin(j*t(i)));
        fpt = fpt - j * j * (harm(j,3)*cos(j*t(i)) + harm(j,4)*sin(j*t(i)));
        gt = gt + j * (harm(j,2)*cos(j*t(i)) - harm(j,1)*sin(j*t(i)));
        gpt = gpt - j * j * (harm(j,1)*cos(j*t(i)) + harm(j,2)*sin(j*t(i)));
    end

    at(i) = ( ( gt * fpt ) - ( ft * gpt ) ) / ( ((ft * ft) + (gt * gt)) ^ (3/2) );
end




figure

[x,y] = efaDraw(harm,nT);

scatter(x,y,5,at,'filled');

axis equal

end