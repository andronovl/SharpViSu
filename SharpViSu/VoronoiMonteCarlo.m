function [ Histograms, intersection ] = VoronoiMonteCarlo( A, BW, iter, signif )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
fov = FOV(A);
p = fov/size(BW,1);% pixel size of the mask in nm
% cutting the eventlist
Anew = parroifilter(A, BW, 1, p);
Voronoi = VorArea_ia( Anew );
%figure; hist(Voronoi{1}(Voronoi{1}<5*median(Voronoi{1}(Voronoi{1}<Inf))),round(2*(length(Voronoi{1}))^(1/3)));
%
Nbins = round(2*(length(Voronoi{1}))^(1/3));
lim = 4*median(Voronoi{1}(Voronoi{1}<Inf));
[counts, centers] = hist(Voronoi{1}(Voronoi{1}<lim),Nbins);

% monte-carlo
if exist ('iter', 'var') && iter > 0
    h = waitbar(0, 'Voronoi Monte Carlo simulation');  
    z = norminv(1-(100-signif)*0.01/2); % z for given significance level (signif %)
    for j = 1:iter
        Ar = zeros(round((fov*fov) * size(Anew,1) / (bwarea(BW) * p^2)), 9);
        Ar(:,4:5) = rand(size(Ar, 1), 2) * fov;
        Anewr = parroifilter(Ar, BW, 1);
        Voronoi_r = VorArea_ia( Anewr );
        [counts_r(j,:), centers] = hist(Voronoi_r{1}(Voronoi_r{1}<lim),Nbins);
        waitbar(j/iter, h, 'Voronoi Monte Carlo simulation');
    end
close(h);

MeanCounts = mean(counts_r);
StdCounts = std(counts_r);

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
[xi, yi] = polyxpoly(centers, counts, centers, MeanCounts);
intersection = [xi, yi];
Histograms = [centers', counts', MeanCounts', (MeanCounts - z * StdCounts)', (MeanCounts + z * StdCounts)'];
else
    Histograms = [centers', counts'];
    intersection = [0, 0];
end

