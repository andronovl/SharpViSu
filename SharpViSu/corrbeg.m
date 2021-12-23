function [Acorr] = corrbeg(A, offset) 
% shifts the events to the beginning of the acquisition using offset as
% the frame-to-frame drift value (for correction of the channel imaged second [488 nm or 532 nm by default])
if ~isempty (A)
    
s = size(offset,1);


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
if s > 2
    Acorr{1} = A(1:middle_event(2),:);
else
    Acorr{1} = A;
end

first_frame = Acorr{1}(1,2);
Acorr{1}(:,4:5) = Acorr{1}(:,4:5) - offset(1,1:2) .* (Acorr{1}(:,2) - first_frame);


sumoffset = zeros(1,2);
if s > 3
for i = 2:s-2
    Acorr{i} = A(middle_event(i)+1:middle_event(i+1), :);
    first_frame = Acorr{i}(1,2);
    Acorr{i}(:,4:5) = Acorr{i}(:,4:5) - 1.5 .* totoffset(1,1:2) - sumoffset(1:2) - offset(i,1:2) .* (Acorr{i}(:,2) - first_frame);
    sumoffset = sumoffset + totoffset(i, :);
end
end

if s > 2
    Acorr{s-1} = A(middle_event(s-1)+1:end, :);
    first_frame = Acorr{s-1}(1,2);
    Acorr{s-1}(:,4:5) = Acorr{s-1}(:,4:5) - 1.5 .* totoffset(1,1:2) - sumoffset(1:2) - offset(s-1,1:2) .* (Acorr{s-1}(:,2) - first_frame);
end

Acorr = cell2mat(Acorr);
else
    Acorr = A;
end