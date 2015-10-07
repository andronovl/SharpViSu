function [ BW, clusters ] = VoronoiSegmentation( AB, thresh, p, limits )
%Performs a segmentation of localization data using Voronoi diagram.
%Voronoi cells with areas smaller than thresh will be kept and merged with
%their neigbours; 
%input: A, localization table in SharpGSDIM format (4 th and 5 th columns = x and y coordinates of single-molecular events); thresh,
%threshold on area; p, pixel size of binary image, limits = [NevMin,
%NevMax] - exclude clusters with less or equal than NevMin or with equal or
%more than NevMax events
%(optional)
%Output: BW, segmented binary image; clusters, cell array with parameters
%of segmented clusters: clusters{:,1} = x and y coordintes of vertices,
%clusters{:,2} = number of localizations in the cluster, clusters{:,3} =
%area of the cluster; , clusters{:,4} = equivalent eccentricity

h = waitbar(0,'Building Voronoi diagram...');

Voronoi = VorArea_ia(AB);
C = Voronoi{3};
V = Voronoi{2};
area = Voronoi{1};
ic = Voronoi{6};
[n, ~] = histc(ic, unique(ic)); % n = number of repeats for each point with index ic.
%Aun = A(Voronoi{5},:); %unique data points, indexes corresponding to the diagram
[B, index] = sortrows(area);
Csort = C(index);
n = n(index);
%Aunsort = Aun(index,:);
Csmat = cat(2, Csort{:});
polygon = [];
surface = [];
repeat = [];

for i = 1:size(Csort, 1)
    polygon1 = i * ones(numel(Csort{i}),1);
    surface1 = B(i) * ones(numel(Csort{i}),1);
    n1 = n(i) * ones(numel(Csort{i}),1);
    polygon = [polygon; polygon1];
    surface = [surface; surface1];
    repeat = [repeat; n1];
    if rem(i,1000) == 0
        waitbar(i/size(Csort, 1))
    end
end
Csmat = [polygon,surface,Csmat',repeat];


%
%thresh = 500; %threshold on area in nm^2
waitbar(0, h, 'Search for neighbors...');

maxi = sum(Csmat(:,2) <= thresh); %index of the biggest cluster to analyse
cluster = cell(Csmat(maxi,1),1);

Csmat1 = Csmat(1:maxi,:);
for i = 1:maxi
    cluster{Csmat(i,1)} = [cluster{Csmat(i,1)}, (Csmat1(Csmat1(:,3) == Csmat(i,3)))']; %in each cell element there are indices of neigbouring poligons (bigger than threshold)
    if rem(i,1000) == 0
        waitbar(i/maxi)
    end
end
UniqueClusters = cellfun(@unique,cluster,'UniformOutput', 0);
close(h);

UniqueClusters2 = VorClusters(UniqueClusters);


% remove excess clusters
h = waitbar(0, 'Removing excess clusters...');

UniqueClusters2 = cellfun(@unique,UniqueClusters2,'UniformOutput', 0);
%NofRepeats = ones(size(UniqueClusters2));
for i = 1:size(UniqueClusters2,1)
    if ~isempty(UniqueClusters2{i})
    for j = i+1:size(UniqueClusters2,1)
        if isequal(UniqueClusters2{j}, UniqueClusters2{i})
            UniqueClusters2{j} = [];
%            NofRepeats(i) = NofRepeats(i) + 1;
        end
    end
    end
    if rem(i,1000) == 0
        waitbar(i/size(UniqueClusters2,1))
    end
end
UniqueClusters2 = UniqueClusters2(~cellfun(@isempty, UniqueClusters2));




% polybool from mapping toolbox
waitbar(0, h, 'Merging neighboring polygons...');

clusters = cell(size(UniqueClusters2, 1),4);
for i = 1:size(UniqueClusters2,1)
    k = find(Csmat(:,1) ==  UniqueClusters2{i}(1))';
    %area1 = Csmat(k(1),2);
    ver = Csmat(k,3);
    Vnew = {{V(ver,1)}, {V(ver,2)}};
    [Vnew{1}, Vnew{2}] = poly2cw(Vnew{1}, Vnew{2});
    Nevents = Csmat(k(1),4);
    if size(UniqueClusters2{i},2) > 1
    for j = 2:size(UniqueClusters2{i},2)
        k = find(Csmat(:,1) ==  UniqueClusters2{i}(j))';
        ver1 = Csmat(k,3);
        Nevents1 = Csmat(k(1),4);
        Nevents = Nevents + Nevents1; % count number of events in cluster
        %area = Csmat(k(1),2);
        Vnew2 = {{V(ver1,1)}, {V(ver1,2)}};
        [Vnew2{1}, Vnew2{2}] = poly2cw(Vnew2{1}, Vnew2{2});
        [Vnew{1}, Vnew{2}] = polybool('union', Vnew{1}, Vnew{2}, Vnew2{1}, Vnew2{2});
        %area1 = area1 + area;
    end
    else
        Vnew{1}{1}(end+1) = Vnew{1}{1}(1);
        Vnew{2}{1}(end+1) = Vnew{2}{1}(1);
    end
    clusters{i,1} = Vnew; % vertices   
    
    %Nevents = round(size(UniqueClusters2{i},2) * area1/realArea);
    clusters{i,2} = Nevents; %number of events
    if rem(i,1000) == 0
        waitbar(i/size(UniqueClusters2,1))
    end
end

% cut clusters based on limit
if exist('limits', 'var')
    clusters = clusters((cell2mat(clusters(:,2)) > limits(1) & cell2mat(clusters(:,2)) < limits(2)),:);
end


% build binary image
%p = 1; %pixelsize in nm

waitbar(0, h, 'Building segmented image...');

s = 18000/p; %size of image
I1 = zeros(s);
for i = 1:size(clusters,1)
    Area = 0;
    SIx = 0;
    SIy = 0;
    SCxkAk = 0;
    SCykAk = 0;
    for j = 1:size(clusters{i,1}{1},2)
        if ispolycw(clusters{i,1}{1}{j}, clusters{i,1}{2}{j}) %reverse vertice order
            [clusters{i,1}{1}{j}, clusters{i,1}{2}{j}] = poly2ccw(clusters{i,1}{1}{j}, clusters{i,1}{2}{j});
        else
            [clusters{i,1}{1}{j}, clusters{i,1}{2}{j}] = poly2cw(clusters{i,1}{1}{j}, clusters{i,1}{2}{j});
        end
        box = [min(clusters{i,1}{1}{j})/p, min(clusters{i,1}{2}{j})/p; max(clusters{i,1}{1}{j})/p, max(clusters{i,1}{2}{j})/p];
        boxint = [floor(box(1,1)-2), floor(box(1,2)-2); ceil(box(2,1)), ceil(box(2,2))];
        mask = poly2mask(clusters{i,1}{1}{j}/p - box(1,1), clusters{i,1}{2}{j}/p - box(1,2), boxint(2,2) - boxint(1,2), boxint(2,1) - boxint(1,1));
        Ak = zeros(size(clusters{i,1}{1}{j},1)-1,1);
        Cxk = zeros(size(clusters{i,1}{1}{j},1)-1,1);
        Cyk = zeros(size(clusters{i,1}{1}{j},1)-1,1);
        Ixk = zeros(size(clusters{i,1}{1}{j},1)-1,1);
        Iyk = zeros(size(clusters{i,1}{1}{j},1)-1,1);
        for k = 1:size(clusters{i,1}{1}{j},1)-1
            Ak(k) = clusters{i,1}{1}{j}(k) * clusters{i,1}{2}{j}(k+1) - clusters{i,1}{1}{j}(k+1) * clusters{i,1}{2}{j}(k);
            Cxk(k) = clusters{i,1}{1}{j}(k) + clusters{i,1}{1}{j}(k+1);
            Cyk(k) = clusters{i,1}{2}{j}(k) + clusters{i,1}{2}{j}(k+1);
            Ixk(k) = clusters{i,1}{2}{j}(k)^2 + clusters{i,1}{2}{j}(k) * clusters{i,1}{2}{j}(k+1) + clusters{i,1}{2}{j}(k+1)^2;
            Iyk(k) = clusters{i,1}{1}{j}(k)^2 + clusters{i,1}{1}{j}(k) * clusters{i,1}{1}{j}(k+1) + clusters{i,1}{1}{j}(k+1)^2;
         end
%         Ak(end) = clusters{i,1}{1}{j}(end) * clusters{i,1}{2}{j}(1) - clusters{i,1}{1}{j}(1) * clusters{i,1}{2}{j}(end);
%         Cxk(end) = clusters{i,1}{1}{j}(end) + clusters{i,1}{1}{j}(1);
%         Cyk(end) = clusters{i,1}{2}{j}(end) + clusters{i,1}{2}{j}(1);
%         Ixk(end) = clusters{i,1}{2}{j}(end)^2 + clusters{i,1}{2}{j}(end) * clusters{i,1}{2}{j}(1) + clusters{i,1}{2}{j}(1)^2;
%         Iyk(end) = clusters{i,1}{1}{j}(end)^2 + clusters{i,1}{1}{j}(end) * clusters{i,1}{1}{j}(1) + clusters{i,1}{1}{j}(1)^2;
        
        A = (1/2) * sum(Ak); %signed area of the (sub)polygon
        Area = Area + A;
        CxkAk = sum( Cxk .* Ak );
        CykAk = sum( Cyk .* Ak );
        SCxkAk = SCxkAk + CxkAk;
        SCykAk = SCykAk + CykAk;
        Ix = (1/12) * sum( Ixk .* Ak );
        Iy = (1/12) * sum( Iyk .* Ak );
        SIx = SIx + Ix;
        SIy = SIy + Iy;
        
        if sign(A) < 0
            mask = -1 * mask;
        end
        I1(boxint(1,2):boxint(2,2)-1, boxint(1,1):boxint(2,1)-1) = I1(boxint(1,2):boxint(2,2)-1, boxint(1,1):boxint(2,1)-1) + mask;

    end
    Cx = (1/(6*Area)) * SCxkAk;
    Cy = (1/(6*Area)) * SCykAk;
    Ixc = SIx - Area * Cy^2;
    Iyc = SIy - Area * Cx^2;
    if Ixc < Iyc
        eccentricity = sqrt(1-Ixc/Iyc);
    else
        eccentricity = sqrt(1-Iyc/Ixc);
    end
    clusters{i,3} = Area; %cluster area
    clusters{i,4} = eccentricity; %cluster equivalent eccentricity
    if rem(i,1000) == 0
        waitbar(i/size(clusters,1))
    end
end
BW = im2bw(I1);

close(h);

end

