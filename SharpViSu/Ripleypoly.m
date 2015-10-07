function [ R, K, L, Lminusr ] = Ripleypoly( A2, BW1, dr, rmax, BW2 )
%Ripley Calculates Ripley's K, L and L-r functions within radii from dr to rmax with
%step dr, for each event (a row) from the table of coordinates A2;
%inside ROI = binary BW image
%BW2 = a ROI from which the A2 was made (optional, BW2 >= BW1)
%p = pixel size in nm,
%input r in nm,
%output in nm

p = 18000/size(BW1, 1); % pixel size in nm

A1 = parroifilter(A2, BW1, 1); % create A1 = BW1(A2)
Nsteps = round(rmax/dr);
dist = pdist2(A1(:,4:5),A2(:,4:5));
K = zeros(Nsteps+1,1);
R = zeros(Nsteps+1,1);
for l = 1:Nsteps+1
    a = dr * (l-1);  
    K(l) = sum(sum(dist <= a, 2) - 1)/size(dist,1);
    R(l) = a;
end
if exist('BW2', 'var')
    area = bwarea(BW2) * p^2;
    dens = size(A2,1)/area;
else
    area = bwarea(BW1) * p^2;
    dens = size(A1,1)/area;
end
K = K/dens;

L = sqrt(K/pi);

Lminusr = L - R;
%Lminusr(:,3) = 1.42*30/(size(A1,1));
%Lminusr(:,4) = -1.42*30/(size(A1,1));
end

