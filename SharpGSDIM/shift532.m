function [ Rnew ] = shift532( R )
%corrects chromatic aberrations by shifting the 532 nm eventlist

% import a 532nm eventlist
if isstruct(R)
R = R.data;
end

Rnew = R;

if isdeployed % Stand-alone mode.
    [~, result] = system('path');
    currentDir = char(regexpi(result, 'Path=(.*?);', 'tokens', 'once'));
else % MATLAB mode.
    currentDir = pwd;
end

if exist([currentDir '\532.mat'], 'file')
    
    tform = importdata([currentDir '\532.mat']);
    Rnew(:,4:5) = 1000 * transformPointsInverse(tform, R(:,4:5)/1000);
    
else
    
green = importdata([currentDir '\532.dat']);
red = importdata([currentDir '\642.dat']);

dif = green - red;
p = [polyfit(green(:,1),dif(:,1),2); polyfit(green(:,2),dif(:,2),2)];

 % shift correction of the data
Rnew(:,4)=R(:,4) - polyval(p(1,:), R(:,4));
Rnew(:,5)=R(:,5) - polyval(p(2,:), R(:,5));

end
end

