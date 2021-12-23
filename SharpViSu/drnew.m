function [I] = drnew(A, p, fov, ko, ke)
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
    fov = [max(A(:,5)), max(A(:,4))];
end
if max(size(fov)) == 1
    fov = [fov, fov];
end
edgesX = 0:p:fov(2);
edgesY = 0:p:fov(1);
I = histcounts2 (A(ko:ke,5), A(ko:ke,4) ,edgesY, edgesX);
end

