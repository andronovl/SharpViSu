function [g, I1] = gimage(A, BW, R, dr, p, ra)


g = g_at_r_ROI( A, BW, R, R + dr/2, ra);

% Keep only useful columns, transform to nm
gnew = g(:,4:5);
gnew(:,3) = g(:,10); % gnew in nm


fov = FOV(A);
pBW = fov/size(BW,1); % pixel size of the mask
% make query matrix M
s = fov/p; % size in pixels
M = zeros(s*s,2); %template for query matrix
a = 1;
for i = 1:s
    for j = 1:s
        M(a,1) = j;
        M(a,2) = i;
        a = a+1;
    end
end
M = M * p; %convert to nanometers

% apply the mask to M
M1 = zeros(size(M,1),9);
M1(:,4:5) = M;
M1 = filterroi( M1, BW, 1, pBW);
M = M1(:,4:5);

% interpolate g(R) to M
if size(gnew,1) > 100
    M(:,3) = griddata(gnew(:,1), gnew(:,2), gnew(:,3), M(:,1), M(:,2),'v4');
else
    M(:,3) = zeros(s*s,1);
end
gim = M;
            
% delete negative values
% gim(gim(:,3)<0, 3) = 0;

% build an image from M(x,y,z)
I1 = zeros(s);
for k = 1:size(gim,1)
    I1(gim(k,2)/p, gim(k,1)/p) = gim(k,3);
end
