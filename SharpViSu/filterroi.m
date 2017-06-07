function [ Anew ] = filterroi( Adata, BW, include, p)
%UNTITLED Summary of this function goes here
% includes only roi or excludes it
% in nm
% BW = mask image with pixel size p
la = size((Adata),1);
s = size(BW, 1);

% find min and max pixels of the mask
[a(:,1), a(:,2)] = find(BW);
minmax(1,:) = min(a); % Xmin Ymin
minmax(2,:) = max(a); % Xmax Ymax
minmax = minmax * p; % in nm;

if include == 1 %keep the roi specified by BW
a = Adata(:,4:5);
Anew = Adata;
z = zeros(1,9);
for k = 1:la
    x = a(k,2);
    y = a(k,1);
    if x >= minmax(1,1) && x <= minmax(2,1) && y >= minmax(1,2) && y <= minmax(2,2)
for m=1:s % x
if (x > p*(m-1))&&(x <= p*m)
for n=1:s % y
    if (y > p*(n-1))&&(y <= p*n)
        if ~BW(m,n)
            Anew(k, :) = z;
        end
        break
    end
end
break
end
end
    else
        Anew(k, :) = z;
    end
end
Anew( ~any(Anew,2), : ) = [];  % delete zero rows
end


if include == 0 %empty the roi specified by BW
a = Adata(:,4:5);
Anew = Adata;
z = zeros(1,9);
for k = 1:la
    x = a(k,2);
    y = a(k,1);
    if x >= minmax(1,1) && x <= minmax(2,1) && y >= minmax(1,2) && y <= minmax(2,2)
for m=1:s % x
if (x > p*(m-1))&&(x <= p*m)
for n=1:s % y
    if (y > p*(n-1))&&(y <= p*n)
        if BW(m,n)
           Anew(k, :) = z;
        end
        break
    end
end
break
end
end
    end
end
Anew( ~any(Anew,2), : ) = [];  % delete zero rows
end
end

