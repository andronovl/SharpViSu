function [ Image ] = slices( A, dz, pixsize, fov, ko, ke )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here.
% in nm
if isstruct(A)
    A = A.data;
end
if ~exist('ko', 'var')
    ko = 1;
end
if ~exist('ke', 'var')
    ke = size(A, 1);
end

if ~exist('fov', 'var') || fov == 0
    fov = FOV(A);
end

A(:,6) = A(:,6) - 500; %shift Z

if max(A(:,6)) - min(A(:,6)) > 1
    A = A(ko:ke, :);
    NZ = round(3000/dz); % number of slices in Z
    Image = zeros(floor(fov/pixsize), floor(fov/pixsize), NZ);

    Anew = cell(NZ,1);
    for i = 1:NZ
        b = i - NZ / 2;
        Anew{i} = A(A(:,6) > dz * (b - 1) & A(:,6) <= dz * b, :);
    end
    %parfor i = 1:NZ
    for i = 1:NZ
        if ~isempty(Anew{i})
            I = drpar(Anew{i}, pixsize, fov);
            Image(:,:,i) = I;
        else
            Image(:,:,i) = zeros(floor(fov/pixsize));
        end
    end
else
    A = A(ko:ke, :);
    Image = drpar(A, pixsize, fov);
end

end