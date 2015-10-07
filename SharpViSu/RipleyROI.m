function [ R, K, L, Lmr ] = RipleyROI( A, BW, dr, rmax, iter, signif )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%   input:
% A = Leica eventlist (.data)
% BW = binary mask = ROI
% dr = step and rmin, nm
% rmax, nm
% iter = number of iterations for Monte Carlo
% signif = significance level in % for critical values

% Output:
% R = vector with increasing radius values, nm
% K(:,1) = Ripley's K-function for corresponding radii from R
% K(:,2) = Lower critical value of Ripley's K-function for corresponding radii from R
% K(:,3) = Upper critical value of Ripley's K-function for corresponding radii from R
% L = Ripley's L-function. columns as for K
% Lmr = L(r) - r. columns as for K
fov = FOV(A);
p = fov/size(BW,1);% pixel size of the mask in nm
% cutting the eventlist
Anew = parroifilter(A, BW, 1, p);

% make new mask
Rerode = round(rmax / p + 1); % radius of the disk for erode
se = strel('disk', Rerode, 4);
BW1 = imerode(BW, se);


% run Ripley for the dataset
% dr = in nm
[ R, K, L, Lmr ] = Ripleypoly( Anew, BW1, dr, rmax, BW );

% run Monte Carlo
if exist ('iter', 'var') && iter > 0
Kr = zeros(round(rmax/dr)+1, iter);
Lr = zeros(round(rmax/dr)+1, iter);
Lmrr = zeros(round(rmax/dr)+1, iter);
z = norminv(1-(100-signif)*0.01/2); % z for given significance level (signif %)
h = waitbar(0, 'Ripley Monte Carlo simulation');        
for k = 1:iter
        Ar = zeros(round((fov*fov) * size(Anew,1) / (bwarea(BW) * p^2)), 9);
        Ar(:,4:5) = rand(size(Ar, 1), 2) * fov;
        Anewr = parroifilter(Ar, BW, 1);
        [ ~, Kr(:,k), Lr(:,k), Lmrr(:,k) ] = Ripleypoly( Anewr, BW1, dr, rmax, BW );
        K(:, 2:3) = [mean(Kr, 2) - z * std(Kr, 0, 2), mean(Kr, 2) + z * std(Kr, 0, 2)];
        L(:, 2:3) = [mean(Lr, 2) - z * std(Lr, 0, 2), mean(Lr, 2) + z * std(Lr, 0, 2)];
        Lmr(:, 2:3) = [mean(Lmrr, 2) - z * std(Lmrr, 0, 2), mean(Lmrr, 2) + z * std(Lmrr, 0, 2)];
        waitbar(k/iter, h, 'Ripley Monte Carlo simulation');
end
        close(h);
end
end

