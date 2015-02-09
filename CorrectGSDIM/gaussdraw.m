function [IN] = gaussdraw(A, p, coeff, fov, ko, ke)
% in nm
if isstruct(A)
Adata = A.data;
else
    Adata = A;
end
if ~exist('coeff', 'var')
    coeff = 240;
end
if ~exist('ko', 'var')
    ko = 1;
end
if ~exist('ke', 'var')
    ke = size(Adata,1);
end

if ~exist('fov', 'var') || fov == 0
    fov = FOV(Adata);
end

s = round(fov / p); % image size in pixels
x = Adata(:,5);
y = Adata(:,4);
if Adata(1,7) ~= 0
    ph = Adata(:,7);
else
    ph = Adata(:,6);
end

% parallel version

step = (ke - ko + 1) / 12;
I = cell(2,1);
xn = cell(2,1);
yn = cell(2,1);
phn = cell(2,1);
for i = 1:12
    kon = round(step * (i-1) + ko);
    ken = round(ke-step * (12-i));
    xn{i} = x(kon:ken);
    yn{i} = y(kon:ken);
    phn{i} = ph(kon:ken);
end
for i = 1:12
    Ia = zeros(s);
    xna = xn{i};
    yna = yn{i};
    pha = phn{i};
    e = size(xna,1);
    for k = 1:e
        for m = 1 : s % x on final image
            if (xna(k) > p * (m-1)) && (xna(k) <= p * m)
                for n = 1 : s % y on final image
                    if (yna(k) > p * (n - 1)) && (yna(k) <= p * n)
                        sigma = coeff/(eps + p * sqrt(pha(k))); % gaussian's sigma in pixels
                        kernel = 2 * ceil((5 * sigma) / 2) - 1; % round kernel to odd number
                        if kernel > 1
                            G = fspecial('gaussian', [kernel kernel], sigma);
                            if ~((m - ( kernel / 2 - 0.5 )) < 1) && ~((m + ( kernel / 2 - 0.5 )) > s) && ~((n - ( kernel / 2 - 0.5 )) < 1) && ~((n + ( kernel / 2 - 0.5 ) ) > s)
                                Ia( m - ( kernel / 2 - 0.5 ) : m + ( kernel / 2 - 0.5 ) , n - ( kernel / 2 - 0.5 ) : n + ( kernel / 2 - 0.5 ) ) = Ia( m - ( kernel / 2 - 0.5 ) : m + ( kernel / 2 - 0.5 ), n - ( kernel / 2 - 0.5 ) : n + ( kernel / 2 - 0.5 ) ) + G;
                            end
                        else
                            Ia(m,n) = Ia(m,n) + 1;
                        end
                        break;
                    end
                end
                break;
            end
        end
    end
    I{i} = Ia;
end

IN = zeros(s);
for i = 1:12
    IN = IN + I{i};
end


% sequential version
%{
IN = zeros(s);
e = size(x,1);
for k = 1:e
    for m = 1 : s % x on final image
        if (x(k) > p * (m-1)) && (x(k) <= p * m)
            for n = 1 : s % y on final image
                if (y(k) > p * (n - 1)) && (y(k) <= p * n)
                    sigma = coeff/(eps + p * sqrt(ph(k))); % gaussian's sigma in pixels
                    kernel = 2 * ceil((5 * sigma) / 2) - 1; % round kernel to odd number
                    if kernel > 1
                        G = fspecial('gaussian', [kernel kernel], sigma);
                        if ~((m - ( kernel / 2 - 0.5 )) < 1) && ~((m + ( kernel / 2 - 0.5 )) > s) && ~((n - ( kernel / 2 - 0.5 )) < 1) && ~((n + ( kernel / 2 - 0.5 ) ) > s)
                        IN( m - ( kernel / 2 - 0.5 ) : m + ( kernel / 2 - 0.5 ) , n - ( kernel / 2 - 0.5 ) : n + ( kernel / 2 - 0.5 ) ) = IN( m - ( kernel / 2 - 0.5 ) : m + ( kernel / 2 - 0.5 ), n - ( kernel / 2 - 0.5 ) : n + ( kernel / 2 - 0.5 ) ) + G;
                        end
                        else
                        IN(m,n) = IN(m,n) + 1;
                    end
                    break;
                end
            end
            break;
        end
    end
end
%}