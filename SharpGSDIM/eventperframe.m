function [ evmax ] = eventperframe( A )
% creates a matrix "evmax" showing the quantity of events in the
% corresponding frame (the frame = (row index - 1))
la = size(A,1);
evmax = zeros(A(end, 2) + 1, 1);
frold = A(1, 2);
for ka = 1:la
fr = A(ka, 2);
if fr ~= frold
    evmax(fr) = A(ka-1, 3);
end
frold = fr;
end
evmax(end) = A(la, 3);
%evperfr = [1:size(evmax,1) evmax];
end

