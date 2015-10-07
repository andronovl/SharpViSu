function [ Anew1 ] = g_at_r_ROI( A, BW, R, dr, ra)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%   input:
% A = eventlist (coordinates in nm)
% BW = binary mask for the ROI
% R = radius for g(R), nm
% dr = size of the disk
% ra = number of random events comparing to the length of A (optional)

% Output:
% A = Ripley's g-function(R) for each event in the 10th column
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
larand = ra * la * arearect/areaBW; % length of the random eventlist
Arand = rand(round(larand),2);
Arand(:,1) = Arand(:,1) .* (minmax(2,2)-minmax(1,2)) + minmax(1,2) - 1;
Arand(:,2) = Arand(:,2) .* (minmax(2,1)-minmax(1,1)) + minmax(1,1) - 1;
Arand1 = zeros(size(Arand, 1),9);
Arand1(:,4:5) = Arand;
Arand1 = parroifilter(Arand1, BW, 1, pBW);
Anew = [Anew; Arand1]; % concatenate the initial and the random eventlists to Anew

end

% make new mask

Rerode = round((R + dr/2) / pBW + 1);
    
BWnew = zeros(size(BW)+2);
BWnew(2:size(BW,1)+1, 2:size(BW,2)+1) = BW; % add 1 empty pixel at borders prior to erosion
se = strel('disk', Rerode, 4);
BW1 = imerode(BWnew, se);
BW1 = BW1(2:end-1, 2:end-1); % cut 1 pixel at borders to retreive the original mask size

% run RipleyK for the dataset
Ri = [R-dr/2, R+dr/2];
Anew1 = parroifilter(Anew, BW1, 1, pBW); % create Anew = BW(A)
dist = pdist2(Anew1(:,4:5),Anew(:,4:5));
K = zeros(size(Anew1,1),2); % K = [K(R-dr/2), K(R+dr/2)] 
for i = 1:2
    K(:,i) = RipleyK( Anew, dist, Ri(i), BW, pBW );
end
gR = 1 / (2 * pi * R * dr) .* (K(:,2) - K(:,1));
Anew1 = [Anew1, gR];
end


function [K] = RipleyK (A, dist, R, BWold, p)
    
K = sum(dist < R, 2) - 1;
area = bwarea(BWold) * p^2;
dens = size(A,1)/area;
K = K/dens;

end


