function [Anew] = parroifilter (Adata, BW, include, p)
% in nm
if ~exist('p', 'var')
    p = FOV(Adata)/size(BW,1);
end
% l = size(Adata, 1);
% if l > 1000
% step = l / 12;
% Ane = cell(12, 1);
% parfor i = 1:12
%     ko = round( ( i - 1 ) * step + 1 );
%     ke = round( i * step );
%     An = Adata( ko:ke , : );
%     Ane {i} = filterroi( An, BW, include, p );
% end
%     kee = 0;
% for i = 1:12
%     ke = size (Ane {i}, 1);
%     if ke > 0
%     Anew ((kee + 1):(kee + ke), :) = Ane {i};
%     end
%     kee = kee + ke; %the last written row
% end
%     Anew = sortrows(Anew);
% else
    Anew = filterroi( Adata, BW, include, p);
% end
end


