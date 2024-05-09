function [td] = EFAThirdDerivative(harm, t)

[nH,~] = size(harm);
nT = length(t);

td = zeros(nT,1);

for i = 1:nT
    ft = 0;
    gt = 0;
    fpt = 0;
    gpt = 0;
    fppt = 0;
    gppt = 0;
    for j = 1:nH
        ft = ft + j * (harm(j,4)*cos(j*t(i)) - harm(j,3)*sin(j*t(i)));
        fpt = fpt - j * j * (harm(j,3)*cos(j*t(i)) + harm(j,4)*sin(j*t(i)));
        fppt = fppt + j * j * j * (harm(j,3)*sin(j*t(i)) - harm(j,4)*cos(j*t(i)));

        gt = gt + j * (harm(j,2)*cos(j*t(i)) - harm(j,1)*sin(j*t(i)));
        gpt = gpt - j * j * (harm(j,1)*cos(j*t(i)) + harm(j,2)*sin(j*t(i)));
        gppt = gppt + j * j * j * (harm(j,1)*sin(j*t(i)) - harm(j,2)*cos(j*t(i)));
    end

    td(i) = ( ( ( gt * fppt ) - ( ft * gppt ) ) / ( ( ( ft * ft ) + ( gt * gt ) ) ^ ( 2 ) ) ) - ( ( 3 * ( ( gt * fpt ) - ( ft * gpt ) ) * ( ( 2 * ft * fpt ) + ( 2 * gt * gpt ) ) ) / ( 2 * ( ( ( ft * ft ) + ( gt * gt ) ) ^ ( 3 ) ) ) );
end




figure

[x,y] = efaDraw(harm,nT);

scatter(x,y,5,td,'filled');

axis equal

end