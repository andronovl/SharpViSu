function [A] = adjustformat(A, format)
% adjusts the format to the standard

if isstruct(A)
A = A.data;
end

if format == 1 % LAS AF
    if sum(A(:, 7)) == 0 % check if old Leica format
        Anew = A;
        Anew(:, 6) = A(:,7);
        Anew(:, 7) = A(:,6);
        A = Anew;
    end
A(:,4:5) = A(:,4:5) * 100; % x, y in nm
A(:,8:9) = A(:,8:9) * 100; % sigmas in nm

elseif format == 2 % QuickPALM 
    p = A(1,5)/A(1,3); % pixel size in nm
    Anew = zeros(size(A,1),9);
    Anew(:,4:6) = A(:,5:7); %X,Y in nm
    Anew(:,7) = A(:,2); %intensity (convert to photons???)
    Anew(:,2) = A(:,15) - 1;%frameID
    Anew(:,8) = sum(A(:,8:9),2) * p; % sigmaX in nm??
    Anew(:,9) = sum(A(:,10:11),2) * p; % sigmaY in nm??
    A = sortrows(Anew);
    fr = A(1,2);
    ev = 1;
    %new eventIDs
    for i = 1:size(A,1)
        if A(i,2) == fr
           A(i,3) = ev;
           ev = ev + 1;
        else
           ev = 1;
           fr = A(i,2);
           A(i,3) = ev;
           ev = ev + 1;
        end
    end
    
elseif format == 3 %RapidSTORM
    Anew = zeros(size(A,1),9);
    Anew(:,4:5) = A(:,1:2); %X,Y in nm
    Anew(:,7) = A(:,4); %intensity (convert to photons???)
    Anew(:,2) = A(:,3); %frameID
    A = sortrows(Anew);
    fr = A(1,2);
    ev = 1;
    %new eventIDs
    for i = 1:size(A,1)
        if A(i,2) == fr
           A(i,3) = ev;
           ev = ev + 1;
        else
           ev = 1;
           fr = A(i,2);
           A(i,3) = ev;
           ev = ev + 1;
        end
    end
elseif format == 4 %Localization microscopy of µManager
    A((A(:,8)>10^6|A(:,8)<10),:) = [];%cut events with >10^6 or <10 photons
    Anew = zeros(size(A,1),9);
    Anew(:,4:5) = A(:,6:7); %X,Y in nm
    Anew(:,7) = A(:,8); %intensity (convert to photons???)
    Anew(:,2) = A(:,3); %frameID
    A = sortrows(Anew);
    fr = A(1,2);
    ev = 1;
    %new eventIDs
    for i = 1:size(A,1)
        if A(i,2) == fr
           A(i,3) = ev;
           ev = ev + 1;
        else
           ev = 1;
           fr = A(i,2);
           A(i,3) = ev;
           ev = ev + 1;
        end
    end
end

% % FALCON
%     Anew = zeros(size(A,1),9);
%     Anew(:,4:5) = A(:,2:3)*100; %X,Y in nm
%     Anew(:,7) = A(:,4); %intensity (convert to photons???)
%     Anew(:,2) = A(:,1); %frameID
%     A = sortrows(Anew);
%     fr = A(1,2);
%     ev = 1;
%     %new eventIDs
%     for i = 1:size(A,1)
%         if A(i,2) == fr
%            A(i,3) = ev;
%            ev = ev + 1;
%         else
%            ev = 1;
%            fr = A(i,2);
%            A(i,3) = ev;
%            ev = ev + 1;
%         end
%     end


end
