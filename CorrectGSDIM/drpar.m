function [I] = drpar(A, p, fov, ko, ke)
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
Anew(:,1) = A(ko:ke,5);
Anew(:,2) = A(ko:ke,4);
edges = 0:p:fov-p;
edges = {edges, edges};
I = hist3 (Anew,'Edges', edges);
end

