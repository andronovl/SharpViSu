function [ B ] = SaveViSP( A, format, FPName )
%Converts data to one of ViSP's formats and saves in in the file FPName
%   Detailed explanation goes here
% cut events with wrong Z
A(A(:,6) > 1000 | A(:,6) < -1000, :) = [];

if format == 1 %.2d
    B = zeros(size(A,1), 4);
    B(:,1:2) = A(:,4:5);
    B(:,3) = A(:,7);
    B(:,4) = A(:,2);
    dlmwrite(FPName, B, '\t');
elseif format == 2 %.2dlp
    B = zeros(size(A,1), 4);
    B(:,1:2) = A(:,4:5);
    B(:,3) = 2.355 * 240 ./ (eps + sqrt(A(:,7))); %FWHM precision
    B(:,4) = B(:,3);
    B(:,5) = A(:,7);
    B(:,6) = A(:,2);
    dlmwrite(FPName, B, '\t');
elseif format == 3 %.3d
    B = zeros(size(A,1), 4);
    B(:,1:2) = A(:,4:5);
    B(:,3) = A(:,6);
    B(:,4) = A(:,7);
    B(:,5) = A(:,2);
    dlmwrite(FPName, B, '\t');
elseif format == 4 %.3dlp
    B = zeros(size(A,1), 4);
    B(:,1:2) = A(:,4:5);
    B(:,3) = A(:,6); %Z
    if A(1,8) ~= 0
        B(:,4) = 2.355 * 2 * A(:,8) ./ (eps + sqrt(A(:,7))); %FWHM precision X
        B(:,5) = 2.355 * 2 * A(:,9) ./ (eps + sqrt(A(:,7))); %FWHM precision Y
    else
        B(:,4) = 2.355 * 240 ./ (eps + sqrt(A(:,7))); %FWHM precision X
        B(:,5) = 2.355 * 240 ./ (eps + sqrt(A(:,7))); %FWHM precision Y
    end
    B(:,6) = 1.5 * (B(:,4) + B(:,5)); %FWHM Z = 3 * (FWHMx + FWHMy) / 2  (needs calibration, depends on Nph and Z)
    B(:,7) = A(:,7);
    B(:,8) = A(:,2);
    dlmwrite(FPName, B, '\t');
end

