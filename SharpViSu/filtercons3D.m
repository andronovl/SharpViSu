function [ Anew1 ] = filtercons3D(A, r, gap, anis, reduce)
%UNTITLED Summary of this function goes here
% in nm
% sums consecutive events within an ellipse, a=rin*sigmax/100, b=rin*sigmay/100 (nm)
% around the first one: 
% new coordinates(x,y,z), sigmas, = weighted mean of the old ones, photon count = sum.
% 
f = waitbar(0, 'Merging...');   

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

if ~exist ('reduce', 'var')
    reduce = false;
end

if length(r) == 1
    r = [r, r, r];
elseif length(r) == 2
    r = [r, max(r)];
end

if reduce == true

la = size((Adata),1);
Anew = Adata;
z = zeros(1,9);
for k = 1:la
    if any(Anew(k,:)) %check if the row is non-zero
    frold = Anew(k,2); %the frame of the first event
    n = 0;
    if anis == 1
        sigma = Anew(k, 8:9)/140; % for sigma from the eventlist
        if sigma(1) * sigma(2) == 0 || sigma(1) * sigma(2) == Inf
            sigma = [1 1];
        end
    else
        sigma = [1 1];
    end

    for i = (k+1):la  
        fr = Adata(i,2);
        if fr <= (frold + 1 + gap) % for gap empty frames
    if (((Anew(i,4)-Anew(k,4))/(r(1)*sigma(1)))^2)+((Anew(i,5)-Anew(k,5))/(r(2)*sigma(2)))^2+(((Anew(i,6)-Anew(k,6))/(r(3)))^2) <= 1
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
    
if mod(k,50000) == 0
waitbar(k/la, f, 'Merging...');
end

end
Anew1 = Anew(any(Anew,2),:);

else  %replace
       
la = size((Adata),1);
Anew = Adata;
z = zeros(1,6);
ind1 = any(Anew,2);
B = [];
for k = 1:la
    if ind1(k) %check if the row is non-zero and not already updated
        n = 0;
        if anis == 1
            sigma = Anew(k, 8:9)/100; % for sigma from the eventlist
            if sigma(1) * sigma(2) == 0
                sigma = [1 1];
            end
        else
                sigma = [1 1];
        end
        ind = k;
        for i = (k+1):la  
            if Adata(i,2) <= (Anew(k,2) + 1 + gap) % for gap empty frames
        if (((Anew(i,4)-Anew(k,4))/(r(1)*sigma(1)))^2)+((Anew(i,5)-Anew(k,5))/(r(2)*sigma(2)))^2+(((Anew(i,6)-Anew(k,6))/(r(3)))^2) <= 1
           B(1,:) = Anew(k,:);
           n = n + 1;
           B(n+1,:) = Anew(i,:);
           ind = [ind, i];
           ind1(i) = false;
        end
            else
                break
            end
        end
        if ~isempty(B)
            Enew = z;
            Enew(7) = sum(B(:,7));
            Enew(4) = sum(B(:,7).*B(:,4)/Enew(7),1);
            Enew(5) = sum(B(:,7).*B(:,5)/Enew(7),1);
            Enew(6) = sum(B(:,7).*B(:,6)/Enew(7),1);
            Enew(8) = sum(B(:,7).*B(:,8)/Enew(7),1);
            Enew(9) = sum(B(:,7).*B(:,9)/Enew(7),1);
            Anew(ind,4:9) = repmat(Enew(:,4:9),length(ind),1);
            sigm = 140/sqrt(sum(B(:,7)));
            Anew(ind,4) = normrnd(Enew(4), sigm, length(ind), 1);
            Anew(ind,5) = normrnd(Enew(5), sigm, length(ind), 1);
            if Enew(6) ~= 0
                Anew(ind,6) = normrnd(Enew(6), 2*sigm, length(ind),1);
            end
            B = [];
        end
    end
    
if mod(k,50000) == 0
waitbar(k/la, f, 'Merging...');
end

end 
Anew1 = Anew;
end

delete(f); %close waitbar

end

