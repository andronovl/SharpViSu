function [Acorr] = corrend(A, offset) 
% shifts the events towards the end of the acquisition using offset as
% the frame-to-frame drift value (for correction of the channel imaged first [642 nm by default])
s = size(offset,1); % number of parts

middle_frame = offset(:,3);
middle_event = zeros(size(middle_frame));
for i = 1:s
[~, middle_event(i)] = min(abs(A(:,2)-middle_frame(i))); %index of the middle event of the eventlist i
end
middle_event = middle_event - 1;

totoffset = zeros(s-1,2);

for i = 1:s-1
offset_length = middle_frame(i+1) - middle_frame(i);
totoffset(i,:) = offset(i,1:2) * offset_length; %total drift between two consecutive eventlists
end

Acorr = cell(s-1,1);

%cut the entire eventlist into correspondings sublists
if s > 2
    Acorr{s-1} = A(middle_event(s-1)+1:end, :); % the last sublist
else
    Acorr{1} = A; % the only sublist in case if s == 2
end

last_frame = Acorr{s-1}(end,2);
Acorr{s-1}(:,4:5) = Acorr{s-1}(:,4:5) + offset(s-1,1:2) .* (last_frame - Acorr{s-1}(:,2));
    
% for n = 1 : ln
% for k = 4 : 5
%     Acorr{s-1}(n,k) = Acorr{s-1}(n,k) + offset(s-1,k-3) * (last_frame - Acorr{s-1}(n,2)); %correction of the last sublist
% end
% end

sumoffset = zeros(1,2);
if s > 3
for i = s-2 : -1 : 2
    Acorr{i} = A(middle_event(i)+1:middle_event(i+1), :);
%     ln = size(Acorr{i},1);
    last_frame = Acorr{i}(end,2);
    Acorr{i}(:,4:5) = Acorr{i}(:,4:5) + 1.5 .* totoffset(s-1,1:2) + sumoffset(1:2) + offset(i,1:2) .* (last_frame - Acorr{i}(:,2));
% for n = 1 : ln
% for k = 4 : 5
%     Acorr{i}(n,k) = Acorr{i}(n,k) + 1.5 * totoffset(s-1,k-3) + sumoffset(k-3) + offset(i,k-3) * (last_frame - Acorr{i}(n,2));
% end
% end
    sumoffset = sumoffset + totoffset(i, :);
end
end

if s > 2
    Acorr{1} = A(1:middle_event(2),:);
%     ln = size(Acorr{1},1);
    last_frame = Acorr{1}(end,2);
    Acorr{1}(:,4:5) = Acorr{1}(:,4:5) + 1.5 .* totoffset(s-1,1:2) + sumoffset(1:2) + offset(1,1:2) .* (last_frame - Acorr{1}(:,2));
% for n = 1 : ln
% for k = 4 : 5
%     Acorr{1}(n,k) = Acorr{1}(n,k) + 1.5 * totoffset(s-1,k-3) + sumoffset(k-3) + offset(1,k-3) * (last_frame - Acorr{1}(n,2));
% end
% end
end

Acorr = cell2mat(Acorr);
