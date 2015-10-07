function [L, I1] = Limage(A, BW, R, p, ra)

L = L_at_r_ROI( A, BW, R, ra);

fov = FOV(A);
pBW = fov/size(BW,1); % pixel size of the mask
% Keep only useful columns, transform to nm
Lnew = L(:,4:5);
Lnew(:,3) = L(:,10); % Lnew in nm
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

% interpolate L(R) to M
if size(Lnew,1) > 100
    M(:,3) = griddata(Lnew(:,1),Lnew(:,2),Lnew(:,3),M(:,1),M(:,2),'v4');
else
    M(:,3) = zeros(s*s,1);
end
Lim = M;
            
% delete negative values
% Lim(Lim(:,3)<0, 3) = 0;

% build an image from M(x,y,z)
I1 = zeros(s);
for k = 1:size(Lim,1)
    I1(Lim(k,2)/p, Lim(k,1)/p) = Lim(k,3);
end

% show images
% figure;
% imagesc(I1);
