function [ I, Avor ] = drawVor1( Area, A1, ic, p, fov )


Areacorr = zeros(size(A1,1),1);
Areacorr(:) = Area(ic,:);
Avor = [A1, Areacorr];
Avor(:,10) = Avor(:,10) .* Avor(:,7) ./ mean(Avor(:,7)); %multiply by number of photons

% drawing
if exist ('p', 'var')
s = round(fov / p); % size in pixels
M = zeros(s * s, 3); %template for query matrix
a = 1;
for i = 1:s
    for j = 1:s
        M(a,1) = j;
        M(a,2) = i;
        a = a+1;
    end
end
M = M * p; %convert to nanometers

% interpolate Anew to M - v4 is too slow, natural OK
M(:,3) = griddata(Avor(:,4), Avor(:,5), Avor(:,10), M(:,1), M(:,2), 'natural');

%delete negative values
M(M(:,3) < 0, 3) = 0;

% build an image from M(x,y,z)
I = zeros(s);
for k = 1 : size(M,1)
    I(M(k,2) / p, M(k,1) / p) = M(k,3);
end
I(isnan(I)) = 0;
else
    I = [];
end
end

