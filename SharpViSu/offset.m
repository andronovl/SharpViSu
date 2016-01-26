function [offset, offset1] = offset(A, corrpixsize, s, method)
% calculates drift by cross-correlation
% input: A - eventlist; corrpixsize - pixel size of the histogram images for
% cross-correlation; s - number of consecutive blocks
% output: offset - frame-to-frame offset, each row corresponds to a pair of consecutive images;
% offset1 - total drift between two consecutive eventlists in nm, each row corresponds to pair of
% consecutive images

l = size(A, 1); % number of rows
border_event = floor(l/s); %the index of the final frame of I1
I12 = cell(s,1);

middle_event = zeros(s,1); %index of the middle event of the eventlist i
middle_frame = zeros(s,1);
for i = 1:s
    ko = 1 + border_event * (i - 1);
    ke = border_event * i;
    if method == 1 %hist
        I12{i} = drpar(A, corrpixsize, 0, ko, ke);
    elseif method == 2 %gauss
        I12{i} = gaussdraw(A, corrpixsize, 240, 0, ko, ke);
    elseif method == 3 %dens
        I12{i} = drawVor(A, corrpixsize, 0, ko, ke);
    end
middle_event(i) = round((ke - ko) / 2 + ko - 1/2);
middle_frame(i) = A(middle_event(i),2);
end

if s > 1
offset = zeros(s-1,2);
offset1 = zeros(s-1,2);
    %parfor i = 1:s-1
    for i = 1:s-1
output = dftregistration(fft2(I12{i+1}),fft2(I12{i}),1000); 
% Citation for this algorithm:
% Manuel Guizar-Sicairos, Samuel T. Thurman, and James R. Fienup, 
% "Efficient subpixel image registration algorithms," Opt. Lett. 33, 
% 156-158 (2008).
corr_offset = [output(4) output(3)];
offset_length = middle_frame(i+1) - middle_frame(i);
frameoffsetpix = corr_offset/offset_length; % frame-to-frame pixel drift
offset1(i,:) = corr_offset * corrpixsize; %total drift between two consecutive eventlists in nm
offset(i,:) = frameoffsetpix * corrpixsize; %frame-to-frame drift in nm
    end
offset1(1,:) = 1.5 .* offset1(1,:);
offset1(end,:) = 1.5 .* offset1(end,:);
else
    offset = NaN;
    offset1 = NaN;
    h = errordlg('Please use at least 2 subsets for drift estimation', 'Incorrect number of subsets');
end
