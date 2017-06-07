function [ fov ] = FOV( A )
% find the size of the field of view with precision of 1 µm

if iscell(A)
    s = max(size(A));
    maxval = zeros(s,1);
    for i = 1:s
        if ~isempty(A{i})
            maxval(i) = max(max(A{i}(:,4:5)));
        else
            maxval(i) = 0;
        end
    end
    maxXY = max(maxval);
else
        if ~isempty(A)
            maxXY = max(max(A(:,4:5)));
        else 
            maxXY = 0;
        end
end
if maxXY > 18000
    fov = ceil(maxXY/1000) * 1000;
else 
    fov = 18000;
end

end

