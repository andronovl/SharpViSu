function [Anew, Bnew] = filtercons3Ddemix(A, B, r, gap, anis)
%UNTITLED Summary of this function goes here
% input = two paired datasets (same number of events in each frame!)
% sums consecutive events within an ellipse, a=rin*sigmax/100, b=rin*sigmay/100 (nm)
% around the first one:
% resuling photon count = sum.
% the other parameters are not refined
%
f = waitbar(0, 'Merging...');

if ~exist ('gap', 'var')
    gap = 0;
end

if ~exist ('anis', 'var')
    anis = 0;
end

if length(r) == 1
    r = [r, r, r];
elseif length(r) == 2
    r = [r, max(r)];
end

la = size((A),1);
Anew = A;
Bnew = B;
ind1 = any(Anew,2);
At = [];
Bt = [];
for k = 1:la
    if ind1(k) %check if the row is non-zero and not already updated
        n = 0;
        if anis == 1
            sigmaA = Anew(k, 8:9)/140; % for sigma from the eventlist
            sigmaB = Bnew(k, 8:9)/140; % for sigma from the eventlist
            if sigmaA(1) * sigmaA(2) == 0
                sigmaA = [1 1];
            end
            if sigmaB(1) * sigmaB(2) == 0
                sigmaA = [1 1];
            end
        else
                sigmaA = [1 1];
                sigmaB = [1 1];
        end
        ind = k;
        for i = (k+1):la  
            if A(i,2) <= (Anew(k,2) + 1 + gap) % for gap empty frames
        if (((Anew(i,4) - Anew(k,4))/(r(1) * sigmaA(1))) ^ 2) + ((Anew(i,5) - Anew(k,5))/(r(2) * sigmaA(2))) ^ 2 + (((Anew(i,6) - Anew(k,6))/(r(3))) ^ 2) <= 1 && (((Bnew(i,4) - Bnew(k,4))/(r(1) * sigmaB(1))) ^ 2)+((Bnew(i,5) - Bnew(k,5))/(r(2) * sigmaB(2))) ^ 2 + (((Bnew(i,6) - Bnew(k,6))/(r(3))) ^ 2) <= 1
           At(1,:) = Anew(k,:);
           Bt(1,:) = Bnew(k,:);
           n = n + 1;
           At(n+1,:) = Anew(i,:);
           Bt(n+1,:) = Bnew(i,:);
           ind = [ind, i];
           ind1(i) = false;
        end
            else
                break
            end
        end
        if ~isempty(At) && ~isempty(Bt)
            Anew(ind,7) = sum(At(:,7));
            Bnew(ind,7) = sum(Bt(:,7));
            At = [];
            Bt = [];
        end
    end


if mod(k,50000) == 0
    waitbar(k/la, f, 'Merging...');
end

end
delete(f); %close waitbar