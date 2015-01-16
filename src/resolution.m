function [ FRC1, resolv ] = resolution(A, n, fig)
%UNTITLED5 Summary of this function goes here
%   Calculates Fourier ring correlation curve 
% n = number of steps (to fov)
% resolv - resolution in nm at FRC = 1/7

h = waitbar(0, 'Calculation of FRC...');
if isstruct(A)
A = A.data;
end
fov = FOV(A);
l = size(A, 1);
Arand = A(randperm(l), :);
lh = round(l/2);
A1 = Arand(1:lh, :);
A2 = Arand(lh+1:end, :);
I1 = drpar(A1, 10, fov);
I2 = drpar(A2, 10, fov);
w = tukeywin(fov/10,0.25);
w2 = w*w';
I11 = I1.*w2;
I12 = I2.*w2;
Fs11 = fftshift(single(fft2(I11)));
Fs12 = fftshift(single(fft2(I12)));
[rr, cc] = meshgrid(1:fov/10);
t = sqrt((rr-fov/20).^2+(cc-fov/20).^2);
t2 = (fov/20) / n;
t3 = n /(fov/20);
FRC = zeros(1, n);
for i = 1:n
C1 = t <= t2 * (i - t3);
C2 = t <= t2 * i;
C3 = C2 - C1;
Fs11c = Fs11 .* C3;
Fs12c = Fs12 .* C3;
%Fs11c( ~any(Fs11c,2), : ) = [];  %delete zero-rows
Fs11c( :, ~any(Fs11c,1) ) = [];  %delete zero-columns
%Fs12c( ~any(Fs12c,2), : ) = [];  %delete zero-rows
Fs12c( :, ~any(Fs12c,1) ) = [];  %delete zero-columns
sumFs1 = sum(sum(Fs11c .* conj(Fs12c)));
Fs11c2 = sqrt(sum(sum(Fs11c .* conj(Fs11c))));
Fs12c2 = sqrt(sum(sum(Fs12c .* conj(Fs12c))));
FRC(i) = sumFs1 / (Fs11c2 .* Fs12c2);
waitbar(i/n);
end

FRC = abs(FRC);
FRC1(:,2) = FRC;
FRC1(:,1) = (1/fov:1/(20*n):1/20);
FRC1(:,3) = smooth(FRC1(:,1),FRC1(:,2),0.35,'loess');
if exist ('fig', 'var')
figure(fig);
plot(FRC1(:,1),FRC1(:,3),'-',FRC1(:,1),FRC1(:,2),'.');
end

for m = 1:n
    if FRC1(m,3) < (1/7)
       a = FRC1(m-2:m+2,:);
       a(:,2) = [];
       break
    end
end
if exist('a', 'var')
p = polyfit(a(:,2),a(:,1),3);
resolv = 1/polyval(p,1/7);
else
    resolv = NaN;
end
close(h);
end

