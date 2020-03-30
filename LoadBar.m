function LoadBar(idx,N,extra)

nDash = 50;

prct = round(idx*100/N);
nPrct = round(idx*nDash/N);

if prct < 10
    S = ['  ',num2str(prct),'% |'];
elseif prct < 100
    S = [' ',num2str(prct),'% |'];
else
    S = [num2str(prct),'% |'];
end

for i = 1:nDash
    if i <= nPrct
        S = strcat(S,'=');
    else
        S = strcat(S,'-');
    end
end
S = strcat(S,'| 100%');
if idx == 0
    n = 0;
else
    n = numel(S) + 1;
end
if nargin > 2
    n = n + extra;
end

fprintf(repmat('\b',1,n))

disp(S);



end

