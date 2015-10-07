function [ UniqueClusters1 ] = VorClusters( UniqueClusters, conver )
%UNTITLED2 Summary of this function goes here
%   recursive function
h = waitbar(0,'Search for the full connectivity...');
if nargin == 1
    conver = false(size(UniqueClusters));
end
%temp = UniqueClusters;
cluster = UniqueClusters;
counter = 0;
for i = (find(~conver))'
    for j = UniqueClusters{i}
        if ~conver(j)
        ind = bsxfun(@eq, cluster{i}, UniqueClusters{j}');
        %if ~all(any(ind,2))
        cluster{i} = [cluster{i}, UniqueClusters{j}(~any(ind,2))];
        %end
        end
    end
    counter = counter + 1;
    if rem(counter,100) == 0
        waitbar(i/(find(~conver, 1, 'last')))
    end
end
%UniqueClusters1 = cellfun(@unique,cluster,'UniformOutput', 0);
conver = cellfun(@isequal, cluster, UniqueClusters);

if any(~conver)
    close(h);
    UniqueClusters1 = VorClusters( cluster, conver );
else
%     for i = 1 : size(UniqueClusters1,1)
%         if ~isempty(UniqueClusters1{i})
%             UniqueClusters1{i} = [i, UniqueClusters1{i}];
%         end
%     end
    UniqueClusters1 = cluster;
close(h);
end
