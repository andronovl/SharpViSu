function [ Histograms, inters ] = VoronoiMonteCarlo( A, BW, iter, signif )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
fov = FOV(A);
p = fov/size(BW,1);% pixel size of the mask in nm
% cutting the eventlist
if size(BW, 2) > 2
    Anew = parroifilter(A, BW, 1, p);
elseif size(BW, 2) == 2
    Anew = A(inROI(images.roi.Freehand(Position = BW), A(:,4), A(:,5)),:);
end
Vor = VorArea_ia( Anew );
S = Vor{1};
%figure; hist(Voronoi{1}(Voronoi{1}<5*median(Voronoi{1}(Voronoi{1}<Inf))),round(2*(length(Voronoi{1}))^(1/3)));
%
% monte-carlo
if exist ('iter', 'var') && iter > 0
    h = waitbar(0, 'Voronoi Monte Carlo simulation');  
    z = norminv(1-(100-signif)*0.01/2); % z for given significance level (signif %)
    for j = 1:iter
        if size(BW, 2) > 2
            Ar = zeros(round((fov*fov) * size(Anew,1) / (bwarea(BW) * p^2)), 9);
            Ar(:,4:5) = rand(size(Ar, 1), 2) * fov;
            Anewr = parroifilter(Ar, BW, 1);
        elseif size(BW, 2) == 2
            Ar = zeros(round((fov*fov) * size(Anew,1) / polyarea(BW(:,1), BW(:,2))), 9);
            Ar(:,4:5) = rand(size(Ar, 1), 2) * fov;
            Anewr = Ar(inROI(images.roi.Freehand(Position = BW), Ar(:,4), Ar(:,5)),:);
        end
         Vorr = VorArea_ia( Anewr );
         Sr{j} = Vorr{1};
        waitbar(j/iter, h, 'Voronoi Monte Carlo simulation');
    end

Nbins = round(2 * (size(S,1)) ^ (1/3));
lim = 3 * median(Sr{1}(Sr{1} < Inf));
[counts, centers] = hist(S(S < lim), Nbins);
counts_r = zeros(iter, Nbins);
for j=1:iter
    [counts_r(j,:), ~] = hist(Sr{j}(Sr{j}<lim),centers);
end
close(h);

MeanCounts = mean(counts_r,1);
StdCounts = std(counts_r,0,1);

%ConfidenceInterval = [MeanCounts', (MeanCounts - z * StdCounts)', (MeanCounts + z * StdCounts)'];

% plot random voronoi
% [vx,vy]=voronoi(ABrand(:,4), ABrand(:,5));
% plot(ABrand(:,4), ABrand(:,5), 'ro', vx, vy, 'b-', 'LineWidth', 1, 'MarkerSize', 6, 'MarkerFaceColor', 'r');
% axis equal
% xlim([0,500])
% ylim([0,500])
% set(gca, 'XTick', [0 100 200 300 400 500], 'YTick', [0 100 200 300 400 500], 'TickLength', [0.02,0.05], 'FontSize', 14, 'YDir', 'reverse');
% 
% print('VorDiagr_rand.tif', '-dtiff', '-r300');
%
% figure; plot(centers, counts, centers, MeanCounts + 1.96*StdCounts, '--m', centers, MeanCounts, '-.g', centers, MeanCounts - 1.96*StdCounts, '--m', 'LineWidth', 1.5);
% xlim([0,4100])
% ylim([0,400])
% set(gca, 'XTick', [0 1000 2000 3000 4000], 'YTick', [0 100 200 300 400], 'TickLength', [0.02,0.05], 'FontSize', 14);
% legend('Data', 'Confidence envelope', 'Mean of random data');
%
%print('curves.tif', '-dtiff', '-r300');

% find intersection
Histograms = [centers', counts', MeanCounts', (MeanCounts - z * StdCounts)', (MeanCounts + z * StdCounts)'];
ind = find((counts - MeanCounts)<0,1);
[xi, yi] = intersection([centers(ind-1), centers(ind)], [counts(ind-1), counts(ind)], [centers(ind-1), centers(ind)], [MeanCounts(ind-1), MeanCounts(ind)]);
inters = [xi, yi];
else
    Histograms = [centers', counts'];
    inters = [0, 0];
end

function [xi, yi] = intersection (x1, y1, x2, y2)
% y = a*x + b;
a1 = (y1(2) - y1(1)) / (x1(2) - x1(1));
a2 = (y2(2) - y2(1)) / (x2(2) - x2(1));
b1 = y1(1) - a1 * x1(1);
b2 = y2(1) - a2 * x2(1);
xi = (b2 - b1) / (a1 - a2);
yi = a1 * xi + b1;
