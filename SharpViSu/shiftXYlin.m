function [shiftx, shifty, shiftz] = shiftXYlin(A, B, frames)

f = waitbar(0, 'Pairing...');
A1n = cell(A(end,2)-A(1,2)+1,1);
B1n = cell(B(end,2)-B(1,2)+1,1);
for i = 1:size(A,1)
    fr = A(i,2);
    if fr >= frames + A(1,2)
        A1n(fr - A(1,2) + 1:end) = [];
        break
    end
    A1n{fr - A(1,2) + 1} = [A1n{fr - A(1,2) + 1}; A(i,:)];
    if mod(i,500000) == 0
        waitbar(i/size(A,1)/3,f);
    end
end
for i = 1:size(B,1)
    fr = B(i,2);
    if fr >= frames + B(1,2)
        B1n(fr - B(1,2) + 1:end) = [];
        break
    end
    B1n{fr - B(1,2) + 1} = [B1n{fr - B(1,2) + 1}; B(i,:)];
    if mod(i,500000) == 0
        waitbar(1/3 + i/size(B,1)/3,f);
    end
end
%
shiftx = cell(size(A1n));
shifty = cell(size(A1n));
shiftz = cell(size(A1n));
for i = 1:size(A1n,1)
    if ~isempty(A1n{i}) && ~isempty(B1n{i})
        [~,I] = pdist2(B1n{i}(:,4:6), A1n{i}(:,4:6), 'euclidean', 'Smallest', 1);
        shiftx{i} = B1n{i}(I, 4) - A1n{i}(:,4);
        shifty{i} = B1n{i}(I, 5) - A1n{i}(:,5);
        shiftz{i} = B1n{i}(I, 6) - A1n{i}(:,6);
    end
    if mod(i,5000) == 0
        waitbar(2/3 + i/size(A1n,1)/3,f);
    end
end
shiftx = cell2mat(shiftx);
shifty = cell2mat(shifty);
shiftz = cell2mat(shiftz);
delete(f);
