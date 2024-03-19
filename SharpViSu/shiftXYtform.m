function [B1new, A1new] = shiftXYtform(A, B, frbeg, frend, maxdist, maxdistZ)

f = waitbar(0, 'Pairing...');

A((A(:,4) < min(B(:,4)) + maxdist) | (A(:,5) < min(B(:,5)) + maxdist) | (A(:,5) > max(B(:,5)) - maxdist) | (A(:,4) > max(B(:,4)) - maxdist), :) = [];

A = A(A(:,2) >= frbeg,:);
B = B(B(:,2) >= frbeg,:);

A1n = cell(A(end,2)-frbeg+1,1);
B1n = cell(B(end,2)-frbeg+1,1);
for i = 1:size(A,1)
    fr = A(i,2);
    if fr >= frend
        A1n(fr - frbeg + 1:end) = [];
        break
    end
    A1n{fr - frbeg + 1} = [A1n{fr - frbeg + 1}; A(i,:)];
    if mod(i,500000) == 0
        waitbar(i/size(A,1)/3,f);
    end
end
for i = 1:size(B,1)
    fr = B(i,2);
    if fr >= frend
        B1n(fr - frbeg + 1:end) = [];
        break
    end
    B1n{fr - frbeg + 1} = [B1n{fr - frbeg + 1}; B(i,:)];
    if mod(i,500000) == 0
        waitbar(1/3 + i/size(B,1)/3,f);
    end
end
%
A1new = cell(size(A1n));
B1new = cell(size(A1n));
for i = 1:size(A1n,1)
    if ~isempty(A1n{i}) && ~isempty(B1n{i})
        [D,I] = pdist2(B1n{i}(:,4:5), A1n{i}(:,4:5), 'euclidean', 'Smallest', 1);
        A1new{i} = A1n{i}(D <= maxdist,:);
        B1new{i} = B1n{i}(I(D <= maxdist),:);    
        if exist("maxdistZ","var")
            del = abs(A1new{i}(:,6) - B1new{i}(:,6)) >  maxdistZ;
            A1new{i}(del, :) = [];
            B1new{i}(del, :) = [];
        end
    end
    if mod(i,5000) == 0
        waitbar(2/3 + i/size(A1n,1)/3,f);
    end
end
A1new = cell2mat(A1new);
B1new = cell2mat(B1new);
delete(f);