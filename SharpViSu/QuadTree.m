function [ I1 ] = QuadTree( A, capacity  )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

Anew(:,1) = A(:,5);
Anew(:,2) = A(:,4);
i = 1;
stop = false(1);
while ~stop
    I{i} = uint8(hist3(Anew, [2^i, 2^i]));
    if max(max(I{i})) <= capacity
        stop = true(1);
    end
    i = i + 1;
end

% resize to the maximum
for i = 1 : numel(I)-1
	I{i} = imresize(I{i}, size(I{numel(I)},1) / size(I{i},1), 'nearest');
end
I = cat(3, I{:});

% delete excess values
for i = 1 : size(I,3)-1
    It = I(:,:,i+1);
    It(I(:,:,i) <= capacity) = 0;
    I(:,:,i+1) = It;
end
clear It;

% delete excess values (>capacity)
I(I>capacity) = 0;

% delete values corresponding to nonzero components of higher levels
for i = size(I,3) : -1 : 2
    Nonzeros = find(I(:,:,i));
    for j = i-1 : -1 : 1
    Ij = I(:,:,j);
    Ij(Nonzeros) = 0;
    I(:,:,j) = Ij;
    end
end

% normalize values on areas
I = single(I);
I = I * (255 / max(max(max(I))));
c = 2;
for i = size(I,3)-1:-1:1
    I(:,:,i) = I(:,:,i) / (2 ^ c);
    c = c + 2;
end

% assemble the image
I1 = sum(I,3);

