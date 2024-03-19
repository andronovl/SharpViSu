function [offset, offset1] = offset3D(A, corrpixsize, s, maxVal, blur)
% calculates drift by cross-correlation
% input: A - eventlist; corrpixsize - pixel size of the histogram images for
% cross-correlation; s - number of consecutive blocks
% output: offset - frame-to-frame offset, each row corresponds to a pair of consecutive images;
% offset1 - total drift between two consecutive eventlists in nm, each row corresponds to pair of
% consecutive images
if corrpixsize(1) > 0 && s > 1 && ~isempty(A)

    if length(corrpixsize) == 1
        corrpixsize(2) = corrpixsize(1);
    end

    if ~exist('maxVal', 'var')
        maxVal = 0;
    end
    if ~exist('blur', 'var')
        blur = 0;
    end

    l = size(A, 1); % number of rows
    border_event = floor(l/s); %the index of the final frame of I1

    middle_event = zeros(s,1); %index of the middle event of the eventlist i
    middle_frame = zeros(s,1);
    for i = 1:s
        ko = 1 + border_event * (i - 1);
        ke = border_event * i;
        I1 = slices(A, corrpixsize(2), corrpixsize(1), 0, ko, ke);
        if maxVal > 0
            I1(I1 > maxVal) = maxVal;
        end
        if ndims(I1) == 3
            if blur > 0 %gauss.blur
                I1 = imgaussfilt3(I1, blur);
            end
            %             Iz = squeeze(sum(I1,[1 2]));
            %             f = fit((1:size(Iz,1))', Iz, 'gauss1');
            %             coeff = coeffvalues(f);
            %             Iz = 0.5 + exp(-(((1:size(Iz,1))-coeff(2))/coeff(3)).^2) ;
            I12xy{i} = squeeze(sum(I1,3));
            I12zx{i} = squeeze(sum(I1(:,round(size(I1,2)/4 - size(I1,2)/10):round(size(I1,2)/4 + size(I1,2)/10),:), 2));
            %             I12zx{i} = I12zx{i}./Iz;
            I12zy{i} = squeeze(sum(I1(round(size(I1,1)/4 - size(I1,1)/10):round(size(I1,1)/4 + size(I1,1)/10),:,:),1));
            %             I12zy{i} = I12zy{i}./Iz;
        else
            if blur > 0 %gauss.blur
                I1 = imgaussfilt(I1, blur);
            end
            I12xy{i} = I1;
        end

        middle_event(i) = round((ke - ko) / 2 + ko - 1/2);
        middle_frame(i) = A(middle_event(i),2);
    end

    if s > 1
        offset = zeros(s,3);
        offset1 = zeros(s-1,3);
        %parfor i = 1:s-1
        for i = 1:s-1
            outputXY = dftregistration_centered(fft2(I12xy{i+1}),fft2(I12xy{i}),100);

            if exist("I12zx","var")
                outputZX = dftregistration_centered(fft2(I12zx{i+1}),fft2(I12zx{i}),100);
                outputZY = dftregistration_centered(fft2(I12zy{i+1}),fft2(I12zy{i}),100);
            else
                outputZX = [0 0 0 0];
                outputZY = [0 0 0 0];
            end
            corr_offset = [outputXY(4) outputXY(3) mean([outputZX(4) outputZY(4)])]; %shift X Y Z
            offset_length = middle_frame(i+1) - middle_frame(i);
            frameoffsetpix = corr_offset/offset_length; % frame-to-frame pixel drift
            offset1(i,:) = corr_offset .* [corrpixsize(1) corrpixsize(1) corrpixsize(2)]; %total drift between two consecutive eventlists in nm
            offset(i,:) = frameoffsetpix .* [corrpixsize(1) corrpixsize(1) corrpixsize(2)]; %frame-to-frame drift in nm
        end
        offset(:,4) = middle_frame;
        offset1(1,:) = 1.5 .* offset1(1,:);
        offset1(end,:) = 1.5 .* offset1(end,:);
    end

else
    offset = NaN;
    offset1 = NaN;
end

% Reference for dftregistration:
% Manuel Guizar-Sicairos, Samuel T. Thurman, and James R. Fienup,
% "Efficient subpixel image registration algorithms," Opt. Lett. 33,
% 156-158 (2008).