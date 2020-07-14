function [rad] = AngleND(V1,V2)

[nR,nC] = size(V1);

if nR > 1 && nC > 1
    error('Input row vectors');
elseif nR > 1
    V1 = V1';
end

[nR,nC] = size(V2);

if nR > 1 && nC > 1
    error('Input row vectors');
elseif nR > 1
    V2 = V2';
end
    


V12 = V1*V2';

V1N = norm(V1);
V2N = norm(V2);

vec_norm = V1N*V2N';
CosTheta = V12./vec_norm;
rad = abs(acos(CosTheta));


end

