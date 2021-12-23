function [shiftxx, shiftyy] = shiftXYpos(A, B, frbeg, frend, maxdist)

f = waitbar(0, 'Pairing...');
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
shiftx = cell(size(A1n));
shifty = cell(size(A1n));
posx = cell(size(A1n));
posy = cell(size(A1n));
for i = 1:size(A1n,1)
    if ~isempty(A1n{i}) && ~isempty(B1n{i})
        [~,Idx] = pdist2(B1n{i}(:,4), A1n{i}(:,4),'euclidean','Smallest',1);
        [~,Idy] = pdist2(B1n{i}(:,5), A1n{i}(:,5),'euclidean','Smallest',1);
        shiftx{i} = B1n{i}(Idx, 4) - A1n{i}(:,4);
        shifty{i} = B1n{i}(Idy, 5) - A1n{i}(:,5);
        posx{i} = A1n{i}(:,4);
        posy{i} = A1n{i}(:,5);
    end
    if mod(i,5000) == 0
        waitbar(2/3 + i/size(A1n,1)/3,f);
    end
end
shiftx = cell2mat(shiftx);
shifty = cell2mat(shifty);
posx = cell2mat(posx);
posy = cell2mat(posy);

shiftxx = sortrows([posx, shiftx]);
shiftyy = sortrows([posy, shifty]);

shiftxx = shiftxx(abs(shiftxx(:,2)) <= maxdist(1),:);
shiftyy = shiftyy(abs(shiftyy(:,2)) <= maxdist(2),:);
delete(f);