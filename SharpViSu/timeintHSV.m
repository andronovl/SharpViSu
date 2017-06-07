function [I1] = timeintHSV(A, p, figurenumber)
% in nm


if isstruct(A)
Adata = A.data;
else
    Adata = A;
end
fov = FOV(Adata);
s = fov / p; % image size in pixels
I = zeros(s,s,2); % blank image
l = length(Adata); % number of rows
for k = 1:l
for m = 1:s % x on final image
if (Adata(k,5) > p * (m-1)) && (Adata(k,5) <= p * m)
for n = 1:s % y on final image
if (Adata(k,4) > p * (n-1)) && (Adata(k,4) <= p * n)
I(m,n,1) = I(m,n,1) + 1;
I(m,n,2) = I(m,n,2) + Adata(k,2);
break
end
end
break
end
end
end


% averaging the time
for m = 1:s
for n = 1:s
if I(m,n,1) ~= 0
I(m,n,2) = I(m,n,2) / I(m,n,1);
end
end
end

% normalization H channel
I(:,:,2) = 2 * I(:,:,2) / (3 * max(max(I(:,:,2))));

% normalization V channel
I(:,:,1) = I(:,:,1) / (max(max(I(:,:,1))));

% combine the layers into a single HSV image
I1 = cat(3, I(:,:,2), ones(s), I(:,:,1));


if exist ('figurenumber', 'var')

figure(figurenumber); imshow(hsv2rgb(I1));

end