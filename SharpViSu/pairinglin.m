function [A1new, B1new] = pairinglin(A, B, prad)

f = waitbar(0, 'Pairing...');
A1n = cell(A(end,2)-A(1,2)+1,1);
B1n = cell(B(end,2)-B(1,2)+1,1);
for i = 1:size(A,1)
    fr = A(i,2);
    A1n{fr - A(1,2) + 1} = [A1n{fr - A(1,2) + 1}; A(i,:)];
    if mod(i,500000) == 0
        waitbar(i/size(A,1)/3,f);
    end
end
for i = 1:size(B,1)
    fr = B(i,2);
    B1n{fr - B(1,2) + 1} = [B1n{fr - B(1,2) + 1}; B(i,:)];
    if mod(i,500000) == 0
        waitbar(1/3 + i/size(B,1)/3,f);
    end
end
%
A1new = cell(A(end,2)-A(1,2)+1,1);
B1new = cell(A(end,2)-A(1,2)+1,1);
for i = 1:min(size(A1n,1),size(B1n,1))
    if ~isempty(A1n{i}) && ~isempty(B1n{i})
    [D,I] = pdist2(B1n{i}(:,4:5), A1n{i}(:,4:5),'euclidean','Smallest',1);
    A1new{i} = A1n{i}(D <=prad,:);
    B1new{i} = B1n{i}(I(D <=prad),:);
    end
    if mod(i,5000) == 0
        waitbar(2/3 + i/size(A1n,1)/3,f);
    end
end
A1new = cell2mat(A1new);
B1new = cell2mat(B1new);
delete(f);