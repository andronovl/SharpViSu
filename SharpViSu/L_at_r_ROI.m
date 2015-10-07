function [ L ] = L_at_r_ROI( A, BW, R, ra)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%   input:
% A = eventlist (coordinates in nm)
% BW = binary mask for the ROI
% R = radius for L(R), nm

% Output:
% L = Ripley's L-function(R) for each event in the 10th column
fov = FOV(A);
pBW = fov/size(BW,1); % pixel size of the mask
% cutting the eventlist
Anew = parroifilter(A, BW, 1, pBW);
if exist ('ra', 'var') && ra ~= 0
% find min and max pixels of the mask
[a(:,1), a(:,2)] = find(BW);
minmax(1,:) = min(a); % Xmin Ymin
minmax(2,:) = max(a); % Xmax Ymax

minmax = minmax * pBW; % in nm;

areaBW = bwarea(BW) * pBW^2; % area of the mask
arearect = (minmax(2,1)-minmax(1,1)) * (minmax(2,2)-minmax(1,2));
la = size(Anew,1); % number of events
larand = ra * la * arearect/areaBW;
Arand = rand(round(larand),2);
Arand(:,1) = Arand(:,1) .* (minmax(2,2)-minmax(1,2)) + minmax(1,2) - 1;
Arand(:,2) = Arand(:,2) .* (minmax(2,1)-minmax(1,1)) + minmax(1,1) - 1;
Arand1 = zeros(size(Arand, 1),9);
Arand1(:,4:5) = Arand;
Anew = [Anew; Arand1];

Anew = parroifilter(Anew, BW, 1, pBW);
end

% make new mask

Rerode = round(R / pBW + 1);
    
BWnew = zeros(size(BW)+2);
BWnew(2:size(BW,1)+1, 2:size(BW,2)+1) = BW;
se = strel('disk', Rerode, 4);
BW1 = imerode(BWnew, se);
BW1 = BW1(2:end-1, 2:end-1);

% run RipleyL for the dataset
L = RipleyL( Anew, BW1, R, BW, pBW );

end

function [Anew] = RipleyL (A, BW, R, BWold, p)
    
Anew = parroifilter(A, BW, 1, p); % create Anew = BW(A2)
dist = pdist2(Anew(:,4:5),A(:,4:5));
K = sum(dist < R, 2) - 1;
area = bwarea(BWold) * p^2;
dens = size(A,1)/area;
K = K/dens;

Anew(:,10) = sqrt(K/pi);

end


