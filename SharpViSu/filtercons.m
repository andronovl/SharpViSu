function [ Anew1 ] = filtercons(A, rin, gap, anis)
%UNTITLED Summary of this function goes here
% in nm
% adds consecutive events within an ellipse, a=rin*sigmax/100, b=rin*sigmay/100 (nm)
% around the first one: 
% new coordinates(x,y,z), sigmas, = mean of old ones, photon count = sum.
% 
if isstruct(A)
Adata = A.data;
else
    Adata = A;
end


if ~exist ('gap', 'var')
    gap = 0;
end


if ~exist ('anis', 'var')
    anis = 0;
end


la = size((Adata),1);
Anew = Adata;
r = rin;
z = zeros(1,9);
for k = 1:la
    if any(Anew(k,:)) %check if the row is non-zero
    frold = Anew(k,2); %the frame of the first event
    n = 0;
    if anis == 1
    sigma = Anew(k, 8:9)/100; % for sigma from the eventlist
    if sigma(1) * sigma(2) == 0
        sigma = [1 1];
    end
    else
        sigma = [1 1];
    end
    for i = (k+1):la  
        fr = Adata(i,2);
        if fr <= (frold + 1 + gap) % for gap empty frames
    if (((Anew(i,4)-Anew(k,4))/(r*sigma(1)))^2)+(((Anew(i,5)-Anew(k,5))/(r*sigma(2)))^2) <= 1
       B(1,:) = Anew(k,:);
       n = n + 1;
       B(n+1,:) = Anew(i,:);
       frold = Anew(i,2);
       Anew(i,:) = z;
    end
        else
            break
        end
    end
    if exist ('B', 'var')
    Enew = z;
    Enew(2:3) = B(1,2:3);
    Enew(4:6) = mean(B(:,4:6),1);
    Enew(7) = sum(B(:,7));
    Enew(8:9) = mean(B(:,8:9),1);
    Anew(k,:) = Enew;
    clear B;
    end
    end
end
Anew1 = Anew(any(Anew,2),:);

end

