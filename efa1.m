function [ coef , x , y, power] = efa1 ( n , ldmk , nHarm , res)
% Function to perform Eliptical Fourier Analysis (EFA) on Landmark data.


%% Perform EFA

% Find the length between each consecutive pair of co-ordinates

t = zeros(n , 1);
t(n,1) = ((ldmk(1,1) - ldmk(n,1))^2 + (ldmk(1,2) - ldmk(n,2))^2)^0.5;

for i = 1:(n - 1)
    
    t(i,1) = ((ldmk(i + 1,1) - ldmk(i,1))^2 + (ldmk(i + 1,2) - ldmk(i,2))^2)^0.5;
    
end

% Calculate A, B, C and D

A = zeros(nHarm,1);
B = zeros(nHarm,1);
C = zeros(nHarm,1);
D = zeros(nHarm,1);
ASum = zeros(nHarm,n);
BSum = zeros(nHarm,n);
CSum = zeros(nHarm,n);
DSum = zeros(nHarm,n);

for j = 1:nHarm
    
    xp = ldmk(1,1) - ldmk(n,1);
    yp = ldmk(1,2) - ldmk(n,2);
   dtp = t(n);
    tp = sum(t(1:n - 1));
   tp1 = sum(t(1:n - 2));
   p2n = 2 * pi * j;
     T = sum(t);
     
     ASum(j,n) = (xp / dtp) * (cos((p2n*tp)/T) - cos((p2n*tp1) / T));
     BSum(j,n) = (xp / dtp) * (sin((p2n*tp)/T) - sin((p2n*tp1) / T));
     CSum(j,n) = (yp / dtp) * (cos((p2n*tp)/T) - cos((p2n*tp1) / T));
     DSum(j,n) = (yp / dtp) * (sin((p2n*tp)/T) - sin((p2n*tp1) / T));
    
    for k = 1:n - 1
       
        xp = ldmk(k + 1,1) - ldmk(k,1);
        yp = ldmk(k + 1,2) - ldmk(k,2);
       dtp = t(k);
        if k == 1
            tp = 0;
            tp1 = sum(t(1:n - 1));
        elseif k == 2
            tp1 = 0;
            tp = t(1);
        elseif k == 3
            tp1 = t(1);
            tp = t(1) + t(2);
        else
            tp1 = sum(t(1:k - 2));
            tp = sum(t(1:k - 1));
        end
       p2n = 2 * pi * j;

         ASum(j,k) = (xp / dtp) * (cos((p2n*tp)/T) - cos((p2n*tp1) / T));
         BSum(j,k) = (xp / dtp) * (sin((p2n*tp)/T) - sin((p2n*tp1) / T));
         CSum(j,k) = (yp / dtp) * (cos((p2n*tp)/T) - cos((p2n*tp1) / T));
         DSum(j,k) = (yp / dtp) * (sin((p2n*tp)/T) - sin((p2n*tp1) / T)); 
        
    end
    
    A(j) = (T / (2 * pi^2 * j^2)) * sum(ASum(j,:));
    B(j) = (T / (2 * pi^2 * j^2)) * sum(BSum(j,:));
    C(j) = (T / (2 * pi^2 * j^2)) * sum(CSum(j,:));
    D(j) = (T / (2 * pi^2 * j^2)) * sum(DSum(j,:));
    
end

the = (2 *(A(1) * B(1) + C(1) * D(1))) / ((A(1)^2) + (C(1)^2) - (B(1)^2) - (D(1)^2));

theta = 0.5 * atan(the);

run = true;

while(run)

    aST = (A(1) * cos(theta)) + (B(1) * sin(theta));
    cST = (C(1) * cos(theta)) + (D(1) * sin(theta));
    
    THE = atan(cST / aST);
    
    E = 1/((aST^2 + cST^2)^0.5);
    
    cT = cos(THE);
    sT = sin(THE);
    
    tMat = [cT sT ; -sT cT];

    mat = zeros(2,2,nHarm);

    for p = 1:nHarm
        
        mat(:,:,p) = (E * (tMat * [A(p) B(p) ; C(p) D(p)] * [cos(p*theta) -sin(p*theta) ; sin(p*theta) cos(p*theta)]));
        
    end

    for q = 1:res

        for Q = 1:nHarm

            X(q,Q) = (mat(1,1,Q) * cos(Q * (q*2*pi)/res)) + (mat(1,2,Q) * sin(Q * (q*2*pi)/res));
            Y(q,Q) = (mat(2,1,Q) * cos(Q * (q*2*pi)/res)) + (mat(2,2,Q) * sin(Q * (q*2*pi)/res));

        end

        x(q) = sum(X(q,:));
        y(q) = sum(Y(q,:));

    end
    
    if max(y) > max(x)
        theta = theta - pi/2;
    else
        run = false;
    end
    
    if mat(1,1,1) < 0
        theta = theta + pi;
        run = true;
    end
    
end


%% Calculate Harmonic Power

power = zeros(nHarm - 1,1);

for s = 2:nHarm
    if s == 2
        power(s - 1) = ((mat(1,1,s)^2 + mat(1,2,s)^2 + mat(2,1,s)^2 + mat(2,2,s)^2) / 2);
    else
        power(s - 1) = ((mat(1,1,s)^2 + mat(1,2,s)^2 + mat(2,1,s)^2 + mat(2,2,s)^2) / 2) + power(s - 2);
    end
end

%% Save Outputs

coef = zeros(nHarm , 4);

for v = 1 : (nHarm)
    
    coef(v,1) = mat(1,1,v);
    coef(v,2) = mat(1,2,v);
    coef(v,3) = mat(2,1,v);
    coef(v,4) = mat(2,2,v);
    
end
end