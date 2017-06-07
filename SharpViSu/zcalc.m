function [A] = zcalc(A, i, m)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
% replaces A(:,6) with Z-values, calculated with polyfit with coefficients p
% m = 0.79 for refr ind = 1.33; m = 1 for refr index = 1.52;
% in nm

if isstruct(A)
A = A.data;
end
if isdeployed % Stand-alone mode.
    [~, result] = system('path');
    currentDir = char(regexpi(result, 'Path=(.*?);', 'tokens', 'once'));
else % MATLAB mode.
    currentDir = pwd;
end
% % 44.6979*2 = axial chromatic aberration
% % 642nm 
% Z10 = 453.0110487344509 - 44.6979 - 0.6127899351919;
% p1 = [-0.932469549777164 -4.391163316595194 -3.101430938103039 -1.296947077620793e+02 4.530110487344509e+02];
% 
% % 532nm 
% Z20 = 315.1101716477435 + 0.6127899351919;
% p2 = [-1.491154933642583 -11.008189410744840 -22.049246098518820 -1.135074622794605e+02 3.151101716477435e+02];
% 
% % 488nm 
% Z30 = 388.2493769662531 + 44.6979 + 45.3045554555715;
% p3 = [1.528435951770679 -6.192066836793174 -3.080485727036640 -1.249186389262468e+02 3.882493769662531e+02];

if i == 1 % 642nm 
    p = importdata([currentDir '\Z642.dat']);
    A(:,6) = m * ( polyval(p, (A(:,8) - A(:,9))) );
end

if i == 2 % 532nm 
    p = importdata([currentDir '\Z532.dat']);
    A(:,6) = m * ( polyval(p, (A(:,8) - A(:,9))) );
end

if i == 3 % 488nm 
    p = importdata([currentDir '\Z488.dat']);
    A(:,6) = m * ( polyval(p, (A(:,8) - A(:,9))) );
end
end


