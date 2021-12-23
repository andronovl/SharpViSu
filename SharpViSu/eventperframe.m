function [ evmax ] = eventperframe( A, smooth )
% creates a matrix "evmax" showing the number of events in the
% corresponding frame (the frame = (row index - 1))
if ~exist('smooth', 'var')
    smooth = 1;
end
if smooth < 1
    smooth = 1;
end
    edges = 0:smooth:max(A(:,2))+1;
    counts = histcounts(A(:,2),edges);
if smooth ~= 1
%     p = polyfit((0:smooth:max(A(:,2)))+0.5, counts, 2);
%     evmax = polyval(p,0.5:(max(A(:,2))+0.5));
    evmax = interp1(1:smooth:length(counts)*smooth, counts, 0.5:(max(A(:,2))+0.5),'linear');
    evmax = evmax ./ smooth;
else
    evmax = counts;
end