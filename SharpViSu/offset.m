function [offset, offset1] = offset(A, corrpixsize, s, max, blur)
% calculates drift by cross-correlation
% input: A - eventlist; corrpixsize - pixel size of the histogram images for
% cross-correlation; s - number of consecutive blocks
% output: offset - frame-to-frame offset, each row corresponds to a pair of consecutive images;
% offset1 - total drift between two consecutive eventlists in nm, each row corresponds to pair of
% consecutive images
if corrpixsize > 0 && s > 1 && ~isempty(A)

if ~exist('max', 'var')
    max = 0;
end
if ~exist('blur', 'var')
    blur = 0;
end
l = size(A, 1); % number of rows
border_event = floor(l/s); %the index of the final frame of I1
I12 = cell(s,1);

middle_event = zeros(s,1); %index of the middle event of the eventlist i
middle_frame = zeros(s,1);
%fov = FOV(A);
fov = 0;
for i = 1:s
    ko = 1 + border_event * (i - 1);
    ke = border_event * i;
        I = drnew(A, corrpixsize, fov, ko, ke);
        if max > 0
        I(I>max) = max;
        end
    if blur > 0 %gauss.blur
        h = fspecial('gaussian',5,blur);
        I = imfilter(I, h, 'replicate','same');
    end
        I12{i} = I;
middle_event(i) = round((ke - ko) / 2 + ko - 1/2);
middle_frame(i) = A(middle_event(i),2);
end

offset = zeros(s,3);
offset1 = zeros(s-1,2);
    %parfor i = 1:s-1
    for i = 1:s-1
output = dftregistration(fft2(I12{i+1}),fft2(I12{i}),100); 
% Citation for this algorithm:
% Manuel Guizar-Sicairos, Samuel T. Thurman, and James R. Fienup, 
% "Efficient subpixel image registration algorithms," Opt. Lett. 33, 
% 156-158 (2008).
corr_offset = [output(4) output(3)];
offset_length = middle_frame(i+1) - middle_frame(i);
frameoffsetpix = corr_offset/offset_length; % frame-to-frame pixel drift
offset1(i,:) = corr_offset * corrpixsize; %total drift between two consecutive eventlists in nm
offset(i,1:2) = frameoffsetpix * corrpixsize; %frame-to-frame drift in nm
    end
offset(:,3) = middle_frame;
offset1(1,:) = 1.5 .* offset1(1,:);
offset1(end,:) = 1.5 .* offset1(end,:);
else
    offset = NaN;
    offset1 = NaN;
    %h = errordlg('Please use at least 2 subsets for drift estimation', 'Incorrect number of subsets');
end
