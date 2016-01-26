function varargout = ClusterViSu(varargin)
% CLUSTERVISU MATLAB code for ClusterViSu.fig
%      CLUSTERVISU, by itself, creates a new CLUSTERVISU or raises the existing
%      singleton*.
%
%      H = CLUSTERVISU returns the handle to a new CLUSTERVISU or the handle to
%      the existing singleton*.
%
%      CLUSTERVISU('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CLUSTERVISU.M with the given input arguments.
%
%      CLUSTERVISU('Property','Value',...) creates a new CLUSTERVISU or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ClusterViSu_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ClusterViSu_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ClusterViSu

% Last Modified by GUIDE v2.5 30-Sep-2015 10:22:09

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ClusterViSu_OpeningFcn, ...
                   'gui_OutputFcn',  @ClusterViSu_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before ClusterViSu is made visible.
function ClusterViSu_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ClusterViSu (see VARARGIN)

% Choose default command line output for ClusterViSu
handles.output = hObject;
if nargin == 5
    AB = varargin{1};
    handles.AB = AB;
    handles.folder = varargin{2};
    BW = cell(2,1);
    handles.BW = BW;
    handles.ABmasked = AB;
    if isempty(AB{1})
        set(handles.popupmenu1, 'Value', 2);
    end
    guidata(hObject, handles);
    updateSRIm(hObject,handles);
%     uiwait(handles.figure1);
else
    handles.AB = cell(2,1);
    handles.ABmasked = cell(2,1);
    guidata(hObject, handles);
end
R = cell(1,2);
K = cell(1,2);
L = cell(1,2);
Lmr = cell(1,2);
Ripley = [R; K; L; Lmr];
handles.Ripley = Ripley;

Histograms = cell(1,2);
intersection = cell(1,2);
VoronoiStats = [ Histograms; intersection ];
handles.VoronoiStats = VoronoiStats;

handles.I0 = cell(2,1);
handles.L = cell(2,2);
handles.Lim = cell(2,3);
handles.LimRGB = cell(2,3);
guidata(hObject, handles);
% Update handles structure

% 

% UIWAIT makes ClusterViSu wait for user response (see UIRESUME)
 %

% --- Outputs from this function are returned to the command line.
function varargout = ClusterViSu_OutputFcn(hObject, eventdata, handles, varargin) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
if nargin == 4
if isfield(handles, 'output')
varargout{1} = handles.output;
end
%varargout{2} = handles.Ripley;
% delete(hObject);
end


% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
updateSRIm(hObject,handles);
updateRipleyGraph(handles);
updateLimRGB(hObject, handles);
thresholding(handles);




function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
updateSRIm(hObject,handles);




function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
updateSRIm(hObject,handles);



% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try
AB = handles.AB;
if ~isempty(AB{1}) || ~isempty(AB{2})
method = get(handles.popupmenu3, 'Value');
if method == 1 % freehand
    h = imfreehand(handles.axes1);
elseif method == 2 % polygon
    h = impoly(handles.axes1);
elseif method == 3 % rectangle
    h = imrect(handles.axes1);
elseif method == 4 % ellipse
    h = imellipse(handles.axes1);
end
% if isequal(get(handles.figure1, 'waitstatus'), 'waiting')
%     % The GUI is still in UIWAIT, use UIRESUME
% temp = 1;
% end
wait(h); 
BW = createMask(h);
delete(h);
channel = get(handles.popupmenu1, 'Value');
    A = AB{channel};
Anew = parroifilter (A, BW, 1);
handles.ABmasked{channel} = Anew;
handles.BW{channel} = BW;
guidata(hObject, handles);
updateSRIm(hObject,handles);
% if exist('temp' , 'var')
% uiwait(handles.figure1);
% end % update the super-resolution image
else
    errordlg('Please load data first!','No data');
end
catch errorObj
    errordlg(getReport(errorObj,'extended','hyperlinks','off'),'Error');
end




% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
channel = get(handles.popupmenu1, 'Value');
AB = handles.AB{channel};
BW = [];
handles.BW{channel} = BW;
handles.ABmasked{channel} = AB;
guidata(hObject, handles);
updateSRIm(hObject,handles);




% --- Executes on button press in pushbutton7.
function pushbutton7_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try
dr = str2double(get(handles.edit3, 'String'));
rmax = str2double(get(handles.edit4, 'String'));
signif = str2double(get(handles.edit5, 'String'));
iter = str2double(get(handles.edit6, 'String'));
channel = get(handles.popupmenu1, 'Value');

AB = handles.AB{channel};
BW = handles.BW{channel};

if get(handles.checkbox3, 'Value')
Ripley = handles.Ripley;
R = Ripley(1,:);
K = Ripley(2,:);
L = Ripley(3,:);
Lmr = Ripley(4,:);
[ R{channel}, K{channel}, L{channel}, Lmr{channel} ] = RipleyROI( AB, BW, dr, rmax, iter, signif );
Ripley = [R; K; L; Lmr];
handles.Ripley = Ripley;
end

if get(handles.checkbox2, 'Value')
    VoronoiStats = handles.VoronoiStats;
    Histograms = VoronoiStats(1,:);
    intersection = VoronoiStats(2,:);
    [ Histograms{channel}, intersection{channel} ] = VoronoiMonteCarlo( AB, BW, iter, signif );
    VoronoiStats = [ Histograms; intersection ];
    handles.VoronoiStats = VoronoiStats;
    set(handles.edit18, 'String', num2str(round(intersection{channel}(1))));
end

guidata(hObject,handles);
updateRipleyGraph(handles);
setVoronoiThreshold(hObject, handles);

catch errorObj
    errordlg(getReport(errorObj,'extended','hyperlinks','off'),'Error');
end


% --- Executes on selection change in popupmenu5.
function popupmenu5_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
updateRipleyGraph(handles);



function edit9_Callback(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
thresholding (handles);



% --- Executes on selection change in popupmenu7.
function popupmenu7_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
updateStats(handles);



% --- Executes on button press in pushbutton11.
function pushbutton11_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
channel = get(handles.popupmenu1, 'Value');
A = handles.AB{channel};
BW = handles.BW{channel};
pI = str2double(get(handles.edit7, 'String'));
noise = str2double(get(handles.edit12, 'String'));
fov = FOV(A);
pBW = fov/size(BW, 1); % pixel size of the mask in nm
R = str2double(get(handles.edit8, 'String'));

I0 = drmasked(handles.ABmasked{channel}, BW, pI, FOV(handles.AB{channel}));
handles.I0{channel} = I0;

mode = get(handles.popupmenu10, 'Value');
if mode == 1 %L(R)
    [L, I1] = Limage(A, BW, R, pI, noise);
elseif mode == 2 % g(R);
    dr = str2double(get(handles.edit3, 'String')); %dr = step for Ripley analysis
    [L, I1] = gimage(A, BW, R, dr, pI, noise);
elseif mode == 3 % Gaussian blur
    I1 = imfilter(I0, fspecial('gaussian', [round(5 * R / pI), round(5 * R / pI)], R / pI), 'same');
elseif mode == 4 % Voronoi diagram
    Voronoi = VorArea(filterroi( A, BW, 1, pBW));
    I1 = Voronoi{1};
elseif mode == 5 % interpolated Voronoi densities
    [ I1, ~ ] = drawVor(filterroi( A, BW, 1, pBW), pI, FOV(handles.AB{channel}));
end
    
    
if mode ~= 3:4 % if not Gaussian blur
% find min and max pixels of the mask
[a(:,1), a(:,2)] = find(BW);
minmax(1,:) = min(a); % Xmin Ymin
minmax(2,:) = max(a); % Xmax Ymax
minmax = minmax * pBW / pI; % in pixels of I1;
I1 = I1(minmax(1,1):minmax(2,1), minmax(1,2):minmax(2,2));
end

if mode ~= 4
    
if exist ('L', 'var')
handles.L{channel, mode} = L;
end
handles.Lim{channel, mode} = I1;
if mode == 5
    mult = str2double(get(handles.edit17, 'String'));
else
    mult = 1;
end
I = I1;
I(I < 0) = 0;
Lmax = max(max(I));
I = uint8(mult * I * 256 / Lmax);
RGB = ind2rgb(I, jet(256));
I0(I0 > 0) = 1;
I0 = ~I0;
RGB = RGB .* cat(3, I0, I0, I0);

handles.LimRGB{channel, mode} = RGB;
else
    handles.LimRGB{channel, mode} = Voronoi;
end
guidata(hObject,handles);
updateLimRGB(hObject, handles);
thresholding (handles);


function updateLimRGB(hObject, handles)
try
channel = get(handles.popupmenu1, 'Value');
mode = get(handles.popupmenu10, 'Value');
RGB = handles.LimRGB{channel, mode};
if ~isempty(RGB)
    if mode ~= 4
        I = handles.Lim{channel, mode};
        Lmax = max(max(I));
        iptsetpref('ImshowAxesVisible','off');
        imshow(RGB, 'Parent', handles.axes3);
        axis(handles.axes3, 'equal');
        set(handles.text52, 'String', Lmax);
    else %tessellations
        cla(handles.axes3);
        C = RGB{3};
        V = RGB{2};
        A = RGB{4};
        hold on
        for i = 1:length(C)
            patch(V(C{i},1), V(C{i},2), 'w', 'Parent', handles.axes3);
        end
            %plot(A(:,4), A(:,5), '.k', 'MarkerSize', 1, 'Parent', handles.axes3);
        hold off
        % find min and max pixels of the mask
        BW = handles.BW{channel};
        [a(:,1), a(:,2)] = find(BW);
        minmax(1,:) = min(a); % Xmin Ymin
        minmax(2,:) = max(a); % Xmax Ymax
        minmax = minmax * FOV(BW) / size(BW, 1); % in nm
        ylim(handles.axes3, [minmax(1,1) minmax(2,1)]);
        xlim(handles.axes3, [minmax(1,2) minmax(2,2)]);
%       xlim(handles.axes3, [0 18000]);
%       ylim(handles.axes3, [0 18000]);
        set(handles.axes3, 'YDir', 'reverse');
    end
else
    cla(handles.axes3);
end

catch errorObj
    errordlg(getReport(errorObj,'extended','hyperlinks','off'),'Error');
end

function updateSRIm(hObject, handles)
% updates figure in the panel "Super-resolution image"
try
channel = get(handles.popupmenu1, 'Value');
A = handles.ABmasked{channel};
BW = handles.BW{channel};
if ~isempty (A)
p = str2double(get(handles.edit1, 'String'));
fov = FOV(handles.AB{channel});
if isempty (BW)
    I = drpar(A, p, fov);
    area = [];
else
    I = drmasked(A, BW, p, fov);
    pBW = fov/size(BW,1);
    area = bwarea(BW) * pBW^2 * 10^-6;
end

Imax = max(max(I));
bright = str2double(get(handles.edit2, 'String'));
I = I * bright * 256 / Imax;
RGB = ind2rgb(uint8(I), hot(256));
iptsetpref('ImshowAxesVisible','on');
imshow(RGB, 'Parent', handles.axes1,'XData', [0 size(I,2) * p * 0.001], 'YData', [0 size(I,1) * p * 0.001]);
set(handles.axes1,'XTickLabel',[], 'XTick',[]);
set(handles.text22, 'String', area);
set(handles.text30, 'String', size(handles.AB{channel}, 1));
set(handles.text33, 'String', size(A, 1));
handles.prevRGB = cell(2,1);
handles.prevRGB{channel} = RGB;
guidata(hObject, handles);
end

catch errorObj
    errordlg(getReport(errorObj,'extended','hyperlinks','off'),'Error');
end



function updateRipleyGraph(handles)
% updates the content of the "Ripley's K analysis" box
try
if isfield(handles, 'Ripley')
    Ripley = handles.Ripley;
    VoronoiStats = handles.VoronoiStats;
    channel = get(handles.popupmenu1, 'Value');
    graph = get(handles.popupmenu5, 'Value'); % the graph to plot
    if ~isempty (Ripley{1,channel})
        R = Ripley{1,channel};
        dr = R(2) - R(1);
        if graph == 1
            K = Ripley{2,channel};
            if size(K, 2) > 1
                plot(handles.axes2, R, K(:,1), R, K(:,2), R, K(:,3));
            else
                plot(handles.axes2, R, K(:,1));
            end
            set(handles.text42, 'String', 'radius, nm');
        elseif graph == 2
            L = Ripley{3,channel};
            if size(L, 2) > 1
                plot(handles.axes2, R, L(:,1), R, L(:,2), R, L(:,3));
            else
                plot(handles.axes2, R, L(:,1));
            end
            set(handles.text42, 'String', 'radius, nm');
        elseif graph == 3
            Lmr = Ripley{4,channel};
            if size(Lmr, 2) > 1
                plot(handles.axes2, R, Lmr(:,1), R, Lmr(:,2), R, Lmr(:,3));
            else
                plot(handles.axes2, R, Lmr(:,1));
            end
            set(handles.text42, 'String', 'radius, nm');
        elseif graph == 4 % g-function = (dk(r)/dr)/(2*pi*r)
            K = Ripley{2,channel};
            eR = R(1:size(R,1)-1)+R(2)/2;
            if size(K, 2) == 3
                g = diff(K)./[(2 * pi * eR * dr), (2 * pi * eR * dr), (2 * pi * eR * dr)];
            elseif size(K, 2) == 1 %no monte carlo
                g = diff(K)./(2 * pi * eR * dr);
            end
            if size(K, 2) > 1
                plot(handles.axes2, eR, g(:,1), eR, g(:,2), eR, g(:,3));
            else
                plot(handles.axes2, eR, g(:,1));
            end
            set(handles.text42, 'String', 'radius, nm');
        end
    end
    if ~isempty (VoronoiStats{1,channel})
            if graph == 5 %Voronoi
            Histograms = VoronoiStats{1,channel};
            if size(Histograms, 2) > 2 % monte carlo
                plot(handles.axes2, Histograms(:,1), Histograms(:,2), Histograms(:,1), Histograms(:,3), Histograms(:,1), Histograms(:,4), Histograms(:,1), Histograms(:,5));
            elseif size(Histograms, 2) == 2 % no monte-carlo
                plot(handles.axes2, Histograms(:,1), Histograms(:,2));
            end
            xlim(handles.axes2, [0, max(Histograms(:,1))]);
            ylim(handles.axes2, [0, max(max(Histograms(:,2:end)))]);
            set(handles.text42, 'String', 'Voronoi polygon area, nm²');
            end
    end
end

catch errorObj
    errordlg(getReport(errorObj,'extended','hyperlinks','off'),'Error');
end

function updateStats(handles)
% updates the histogram for ClusterViSu' statistics
try
mode = get(handles.popupmenu10, 'Value');
channel = get(handles.popupmenu1, 'Value');
if isfield(handles, 'stats')
p = str2double(get(handles.edit7, 'String')); % pixel size in nm
stats = handles.stats{channel};
show = get(handles.popupmenu7, 'Value');
bins = str2double(get(handles.edit25, 'String'));
if show == 1 % cluster area
    if mode == 4
        area = stats(:,2); % cluster area in nm^2
    else
        area = cat(1, stats.Area) * p^2; % cluster area in nm^2
    end
    hist(area, bins, 'Parent', handles.axes5);
    set(handles.text26, 'String', mean(area));
    set(handles.text16, 'String', median(area));
    set(handles.text18, 'String', std(area));
    set(handles.text46, 'String', 'cluster area, nm²');
elseif show == 2 % diameter
    if mode == 4
        diam = 2 * sqrt(stats(:,2)/pi);
    else
        diam = cat(1, stats.EquivDiameter) * p; % cluster diameter in nm
    end
    hist(diam, bins, 'Parent', handles.axes5);
    set(handles.text26, 'String', mean(diam));
    set(handles.text16, 'String', median(diam));
    set(handles.text18, 'String', std(diam));
    set(handles.text46, 'String', 'cluster diameter, nm');
elseif show == 3 % eccentricity
    if mode == 4
        eccentr = stats(:,3);
    else
        eccentr = cat(1, stats.Eccentricity); % cluster area in nm^2
    end
    hist(eccentr, bins, 'Parent', handles.axes5);
    set(handles.text26, 'String', mean(eccentr));
    set(handles.text16, 'String', median(eccentr));
    set(handles.text18, 'String', std(eccentr));
    set(handles.text46, 'String', 'cluster eccentricity');
elseif show == 4 % events per cluster
    if mode == 4
        NBevents = stats(:,1);
    else
        I0 = handles.I0{channel};
        s = size(stats, 1);
        NBevents = zeros(s,1);
        for i = 1:s
            ind = stats(i).PixelIdxList;
            NBevents(i) = sum(I0(ind));
        end
    end
    hist(NBevents, bins, 'Parent', handles.axes5);
    set(handles.text26, 'String', mean(NBevents));
    set(handles.text16, 'String', median(NBevents));
    set(handles.text18, 'String', std(NBevents));
    set(handles.text46, 'String', 'number of events in cluster');
elseif show == 5 % density
    if mode == 4
        Dens = 1000 * stats(:,1) ./ stats(:,2); % events/µm^2
    else
        I0 = handles.I0{channel};
        s = size(stats, 1);
        Dens = zeros(s,1);
        for i = 1:s
            ind = stats(i).PixelIdxList;
            NBeventsi = sum(I0(ind));
            area = stats(i).Area;
            Dens(i) = NBeventsi/(area * (p * 0.001)^2); % events/µm^2
        end
    end
    hist(Dens, bins, 'Parent', handles.axes5);
    set(handles.text26, 'String', mean(Dens));
    set(handles.text16, 'String', median(Dens));
    set(handles.text18, 'String', std(Dens));
    set(handles.text46, 'String', 'density of events in cluster, 1/µm');
end
end

catch errorObj
    errordlg(getReport(errorObj,'extended','hyperlinks','off'),'Error');
end


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isequal(get(hObject, 'waitstatus'), 'waiting')
    % The GUI is still in UIWAIT, use UIRESUME
    uiresume(hObject);
else
    % The GUI is no longer waiting, just close it
    delete(hObject);
end


function [I] = drmasked (A, BW, p, fov)
% draws the histogram image only in the rectangle circumscribing BW
Anew(:,1) = A(:,5);
Anew(:,2) = A(:,4);
pBW = fov/size(BW,1); % pixel size of the mask in nm
% find min and max pixels of the mask
[a(:,1), a(:,2)] = find(BW);
minmax(1,:) = min(a); % Xmin Ymin
minmax(2,:) = max(a); % Xmax Ymax
minmax = minmax * pBW; % in nm;

edgesx = minmax(1,1)-1:p:minmax(2,1);
edgesy = minmax(1,2)-1:p:minmax(2,2);
edges = {edgesx, edgesy};
I = hist3 (Anew,'Edges', edges);


function thresholding (handles)
try
    
thresh = str2double(get(handles.edit9, 'String')); %threshold
minima = str2double(get(handles.edit13, 'String')); %minimum for h-minima transform
water = get(handles.checkbox1, 'Value'); % watershed
excl = str2double(get(handles.edit15, 'String')); %exclude clusters smaller than excl pixels
exclmax = str2double(get(handles.edit16, 'String')); %exclude clusters bigger than excl pixels
limits = [excl, exclmax];
channel = get(handles.popupmenu1, 'Value');
mode = get(handles.popupmenu10, 'Value');
if mode == 5
    thresh = 1/thresh;
end
if mode ~= 4
    Binary = handles.Lim{channel, mode};
    if ~isempty (Binary)
        Binary(Binary < thresh) = 0;
        Binary(Binary >= thresh) = 1;
        Binary = im2bw(Binary, 0.5);
    if water
        D = -bwdist(~Binary);
        D = imhmin(D, minima);
        D(~Binary) = -Inf;
        L = watershed(D);
        Binary(L == 0) = 0;
    end
    CC = bwconncomp(Binary);
    %CCraw = CC;
    % delete ClusterViSu at borders
    % delete zero rows and columns form the mask
    BW = handles.BW{channel};
    BW( ~any(BW,2), : ) = [];  %rows
    BW( :, ~any(BW,1) ) = [];  %columns
    Isize = CC.ImageSize;
    BW = imresize(BW, Isize);
    BWnew = zeros(size(BW)+2);
    BWnew(2:size(BW,1)+1, 2:size(BW,2)+1) = BW;
    se = strel('square', 3);
    BWer = imerode(BWnew, se);
    BWer = BWer(2:end-1, 2:end-1);

    obj = CC.PixelIdxList;
    c = 0;
    for i = 1:size(obj, 2)
        for k = 1:size(obj{i-c})
        if size(obj{i-c},1) < excl+1 || size(obj{i-c},1) > exclmax-1 || BWer(obj{i-c}(k)) == 0 % delete regions containing only one pixel || delete regions touching the border of the mask
        % if BWer(obj{i-c}(k)) == 0 % delete regions touching the border of the mask
            obj(i-c) = [];
            c = c+1;
            break
        end
        end        
    end
    CC.PixelIdxList = obj;
    CC.NumObjects = size(obj,2);

    Binary = zeros(size(Binary));
    Binary(cat(1, obj{:})) = 1;
    Binary = im2bw(Binary, 0.5);
    end
else
    A = handles.AB{channel};
    BW = handles.BW{channel};
    fov = FOV(A);
    pBW = fov/size(BW, 1);
    pI = str2double(get(handles.edit7, 'String'));
    [ Binary, clusters ] = VoronoiSegmentation( filterroi( A, BW, 1, pBW), thresh, pI, limits );
    [a(:,1), a(:,2)] = find(BW);
    minmax(1,:) = min(a); % Xmin Ymin
    minmax(2,:) = max(a); % Xmax Ymax
    minmax = minmax * pBW / pI; % in pixels of I1;
    Binary = Binary(minmax(1,1)-3:minmax(2,1)-3, minmax(1,2)-3:minmax(2,2)-3);
    ClustersMat = [clusters{:,2}; clusters{:,3}; clusters{:,4}]'; 
    CC = bwconncomp(Binary);
end

if ~isempty (Binary)
I0 = handles.I0{channel};
I = I0;
I0(I0 > 0) = 1;
I0 = im2bw(I0, 0.5);
RGB = hsv2rgb(cat(3, ones(size(I0)), I0, I0 | Binary));
iptsetpref('ImshowAxesVisible','off');
imshow(RGB, 'Parent', handles.axes4);

STATS = regionprops(CC, 'All');
if mode == 4
    clusterNb = size(ClustersMat, 1);
    NBevents = sum(ClustersMat(:,1));
    Area = sum(ClustersMat(:,2));
    set(handles.text44, 'String', 0.1 * round(10^-3 * Area / str2double(get(handles.text22, 'String'))));
else
clusterNb = size(STATS,1);
% for clustered events
%STATSraw = regionprops(CCraw, 'Area', 'PixelIdxList');
NBevents = sum(I(cat(1, STATS.PixelIdxList)));
Area = sum(cat(1, STATS.Area));
%Density = NBevents/Area;
%clusterrawNb = size(STATSraw,1);
p = str2double(get(handles.edit7, 'String'));
set(handles.text44, 'String', 0.1 * round(1000 * Area * (p * 0.001)^2/str2double(get(handles.text22, 'String'))));
end
set(handles.text21, 'String', clusterNb);
set(handles.text54, 'String', clusterNb/str2double(get(handles.text22, 'String')));
set(handles.text36, 'String', NBevents);
set(handles.text38, 'String', 0.1 * round(1000 * NBevents/str2double(get(handles.text33, 'String'))));

handles.Binary = Binary;
handles.segmRGB = cell(2,1);
handles.segmRGB{channel} = RGB;
handles.stats = cell(2,1);
if mode == 4
    handles.stats{channel} = ClustersMat;
else
    handles.stats{channel} = STATS;
end
guidata(handles.figure1, handles);
updateStats(handles);
end

catch errorObj
    errordlg(getReport(errorObj,'extended','hyperlinks','off'),'Error');
end



function edit25_Callback(hObject, eventdata, handles)
% hObject    handle to edit25 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
updateStats(handles);



% --- Executes on button press in checkbox1.
function checkbox1_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
thresholding (handles);



function edit13_Callback(hObject, eventdata, handles)
% hObject    handle to edit13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
thresholding (handles);



function edit15_Callback(hObject, eventdata, handles)
% hObject    handle to edit15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
thresholding (handles);



function edit16_Callback(hObject, eventdata, handles)
% hObject    handle to edit16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
thresholding (handles);


% --------------------------------------------------------------------
function Load_mask_Callback(hObject, eventdata, handles)
% hObject    handle to Load_mask (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
channel = get(handles.popupmenu1, 'Value');
text{1} = '642-nm mask';
text{2} = '488/532-nm mask';
if isfield(handles, 'folder')
    pathname = handles.folder;
    [filename, pathname] = uigetfile({'*.tif'; '*.png'; '*.bmp'}, text{channel}, pathname);
else
[filename, pathname] = uigetfile({'*.tif'; '*.png'; '*.bmp'}, text{channel});
end
handles.BW = cell(2,1);
if pathname ~= 0
handles.BW{channel} = imread([pathname filename]);
handles.folder = pathname;
A = handles.AB{channel};
Anew = parroifilter (A, handles.BW{channel}, 1);
handles.ABmasked{channel} = Anew;
end

%clear old data
guidata(hObject,handles);
updateSRIm(hObject, handles);

% --------------------------------------------------------------------
function xytable_Callback(hObject, eventdata, handles)
% hObject    handle to xytable (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
channel = get(handles.popupmenu1, 'Value');
text{1} = '642-nm eventlist';
text{2} = '488/532-nm eventlist';
if isfield(handles, 'folder')
    pathname = handles.folder;
    [filename, pathname] = uigetfile({'*.ascii'; '*.txt'; '*.dat'}, text{channel}, pathname);
else
[filename, pathname] = uigetfile({'*.ascii'; '*.txt'; '*.dat'}, text{channel});
end
if pathname ~= 0
A = importdata([pathname filename]);

if isstruct(A)
A = A.data;
end

% adjust the format to the standard
Ast = zeros(size(A,1), 9);
Ast(:,4:5) = A;
A = Ast;

handles.folder = pathname;
handles.AB{channel} = A;
handles.ABmasked{channel} = A;
end
handles.BW = cell(2,1);

cleardata(handles);
guidata(hObject,handles);
updateSRIm(hObject,handles);


% --------------------------------------------------------------------
function Save_Callback(hObject, eventdata, handles)
% hObject    handle to Save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try
pathname = handles.folder;
dir = uigetdir(pathname, 'Folder to save the data');
dir = [dir '\'];
channel = get(handles.popupmenu1, 'Value');

% save the preview
prevRGB = handles.prevRGB{channel};
pixsize = get(handles.edit1, 'String');
bright = get(handles.edit2, 'String');
if channel == 1 %red eventlist
FPName=[dir 'prev_red_' pixsize '_' bright '.tif'];
elseif channel == 2 %green eventlist
FPName=[dir 'prev_green_' pixsize '_' bright '.tif'];
end
imwrite(prevRGB, FPName, 'Compression', 'lzw');

% save the cluster map image
if get(handles.popupmenu10, 'Value') == 1 % L(R)
    base = 1;
    basestr = 'L(R)';
elseif get(handles.popupmenu10, 'Value') == 2 % g(R)
    base = 2;
    basestr = 'g(R)';
elseif get(handles.popupmenu10, 'Value') == 3 % Gaussian blur
    base = 3;
    basestr = 'G';
elseif get(handles.popupmenu10, 'Value') == 4 % Voronoi diagram
    base = 4;
    basestr = 'VorDiagram';
elseif get(handles.popupmenu10, 'Value') == 5 % Voronoi density map
    base = 5;
    basestr = 'VorDensity';
end
if ~isempty(handles.LimRGB{channel, base})
Lim = handles.LimRGB{channel, base};
bright = get(handles.edit10, 'String');
pixsize = get(handles.edit7, 'String');
rad = get(handles.edit8, 'String');
noise = get(handles.edit12, 'String');
maxi = get(handles.text52, 'String');
if channel == 1 %red eventlist
FPName=[dir basestr '_red_' bright '_' pixsize '_' rad '_' noise '_' maxi '.tif'];
elseif channel == 2 %green eventlist
FPName=[dir basestr '_green_' bright '_' pixsize '_' rad '_' noise '_' maxi '.tif'];
end
if base ~=4
    imwrite(Lim, FPName, 'Compression', 'lzw');
else
RGB = handles.LimRGB{channel, base};
if ~isempty(RGB)
    %axis(handles.axes3, 'equal');
    h3 = figure('visible', 'off');
    C = RGB{3};
    V = RGB{2};
    A = RGB{4};
    hold on
    for i = 1:length(C)
        ppatch = patch(V(C{i},1), V(C{i},2), 'w');
    end
        plot(A(:,4), A(:,5), '.b', 'MarkerSize', 2);
    hold off
    set(ppatch, 'MarkerSize', 3);
    ax = get(ppatch, 'Parent');
    axis(ax, 'equal');
    % find min and max pixels of the mask
    BW = handles.BW{channel};
    [a(:,1), a(:,2)] = find(BW);
    minmax = zeros(2);
    minmax(1,:) = min(a); % Xmin Ymin
    minmax(2,:) = max(a); % Xmax Ymax
    minmax = minmax * FOV(BW) / size(BW, 1); % in nm
    ylim([minmax(1,1) minmax(2,1)]);
    xlim([minmax(1,2) minmax(2,2)]);
    set(ax, 'YDir', 'reverse');
    if channel == 1 %red eventlist
    FPName=[dir basestr '_red'];
    elseif channel == 2 %green eventlist
    FPName=[dir basestr '_green'];
    end
        % set style for graphs
myStyle = hgexport('factorystyle');
myStyle.Format = 'png';
myStyle.Width = 15;
myStyle.Height = 15;
myStyle.Resolution = 300;
myStyle.Units = 'inch';
myStyle.FixedFontSize = 4;

    hgexport(h3, FPName, myStyle, 'Format', 'png');
    close(h3);
end
end
end

% save the segmented image
if isfield(handles, 'segmRGB') && ~isempty(handles.segmRGB{channel})
segmRGB = handles.segmRGB{channel};
pixsize = get(handles.edit7, 'String');
rad = get(handles.edit8, 'String');
noise = get(handles.edit12, 'String');
thresh = get(handles.edit9, 'String');
watersh = get(handles.checkbox1, 'Value');
if watersh
    hmin = get(handles.edit13, 'String');
    waterstr = ['_w' hmin];
else
    waterstr = [];
end
exclmin = get(handles.edit15, 'String');
exclmax = get(handles.edit16, 'String');
if str2double(exclmin) ~=0 || str2double(exclmax) < Inf
    strexcl = ['_excl' exclmin '-' exclmax];
else
    strexcl = [];
end
if channel == 1 %red eventlist
FPName=[dir 'segm_' basestr '_red_' pixsize '_' rad '_' noise '_' thresh waterstr strexcl '.tif'];
elseif channel == 2 %green eventlist
FPName=[dir 'segm_' basestr '_green_' pixsize '_' rad '_' noise '_' thresh waterstr strexcl '.tif'];
end
imwrite(segmRGB, FPName, 'Compression', 'lzw');
end

% save Ripley functions
if isfield(handles, 'Ripley') && ~isempty(handles.Ripley{1, channel})
    
    % set style for graphs
myStyle = hgexport('factorystyle');
myStyle.Format = 'png';
myStyle.Width = 6.41;
myStyle.Height = 4;
myStyle.Resolution = 300;
myStyle.Units = 'inch';
myStyle.FixedFontSize = 8;

% get data
step = get(handles.edit3, 'String');
rmax = get(handles.edit4, 'String');
conf = get(handles.edit5, 'String');
iter = get(handles.edit6, 'String');

    % save K graph    
if channel == 1 %red eventlist
FPName=[dir 'K_red_' step '_' rmax '_' conf '_' iter '.png'];
elseif channel == 2 %green eventlist
FPName=[dir 'K_green_' step '_' rmax '_' conf '_' iter '.png'];
end
h3 = figure('visible', 'off');
Ripley = handles.Ripley;
R = Ripley{1,channel};
K = Ripley{2,channel};
if size(K, 2) > 1
    plot(R, K(:,1), R, K(:,2), R, K(:,3));
else
    plot(R, K(:,1));
end
xlabel('r(nm)')
ylabel('K(r)')
hgexport(h3, FPName, myStyle, 'Format', 'png');
close(h3);

% save L graph  
if channel == 1 %red eventlist
FPName=[dir 'L_red_' step '_' rmax '_' conf '_' iter '.png'];
elseif channel == 2 %green eventlist
FPName=[dir 'L_green_' step '_' rmax '_' conf '_' iter '.png'];
end
h3 = figure('visible', 'off');
L = Ripley{3,channel};
if size(L, 2) > 1
    plot(R, L(:,1), R, L(:,2), R, L(:,3));
else
    plot(R, L(:,1));
end
xlabel('r(nm)')
ylabel('L(r)')
hgexport(h3, FPName, myStyle, 'Format', 'png');
close(h3);

% save Lmr graph  
if channel == 1 %red eventlist
FPName=[dir 'Lmr_red_' step '_' rmax '_' conf '_' iter '.png'];
elseif channel == 2 %green eventlist
FPName=[dir 'Lmr_green_' step '_' rmax '_' conf '_' iter '.png'];
end
h3 = figure('visible', 'off');
Lmr = Ripley{4,channel};
if size(Lmr, 2) > 1
    plot(R, Lmr(:,1), R, Lmr(:,2), R, Lmr(:,3));
else
    plot(R, Lmr(:,1));
end
xlabel('r(nm)')
ylabel('L(r) - r')
hgexport(h3, FPName, myStyle, 'Format', 'png');
close(h3);

% save g(r) graph  
if channel == 1 %red eventlist
FPName=[dir 'g_red_' step '_' rmax '_' conf '_' iter '.png'];
elseif channel == 2 %green eventlist
FPName=[dir 'g_green_' step '_' rmax '_' conf '_' iter '.png'];
end
h3 = figure('visible', 'off');
K = Ripley{2,channel};
dr = R(2) - R(1);
eR = R(1:size(R,1)-1)+R(2)/2;

if size(K, 2) > 1
    g = diff(K)./[(2 * pi * eR * dr), (2 * pi * eR * dr), (2 * pi * eR * dr)];
    plot(eR, g(:,1), eR, g(:,2), eR, g(:,3));
else
    g = diff(K)./(2 * pi * eR * dr);
    plot(eR, g(:,1));
end
xlabel('r(nm)')
ylabel('g(r)')
hgexport(h3, FPName, myStyle, 'Format', 'png');
close(h3);

end


% save Voronoi cell functions
if isfield(handles, 'VoronoiStats') && ~isempty(handles.VoronoiStats{1, channel})
    
% set style for graphs
myStyle = hgexport('factorystyle');
myStyle.Format = 'png';
myStyle.Width = 6.41;
myStyle.Height = 4;
myStyle.Resolution = 300;
myStyle.Units = 'inch';
myStyle.FixedFontSize = 8;

% get data
conf = get(handles.edit5, 'String');
iter = get(handles.edit6, 'String');
VoronoiStats = handles.VoronoiStats;
Histograms = VoronoiStats{1,channel};
intersection = VoronoiStats(2,:);
inters = num2str(round(intersection{channel}(1)));

% save cell areas histogram    
if channel == 1 %red eventlist
    FPName=[dir 'Vor_red_' conf '_' iter '_' inters '.png'];
elseif channel == 2 %green eventlist
    FPName=[dir 'Vor_green_' conf '_' iter '_' inters '.png'];
end
h3 = figure('visible', 'off');
if size(Histograms, 2) > 2 % monte carlo
    plot(handles.axes2, Histograms(:,1), Histograms(:,2), Histograms(:,1), Histograms(:,3), Histograms(:,1), Histograms(:,4), Histograms(:,1), Histograms(:,5));
elseif size(Histograms, 2) == 2 % no monte-carlo
    plot(handles.axes2, Histograms(:,1), Histograms(:,2));
end
xlim(handles.axes2, [0, max(Histograms(:,1))]);
ylim(handles.axes2, [0, max(max(Histograms(:,2:end)))]);
xlabel('Voronoi polygon area (nm^2)')
ylabel('Occurrence')
hgexport(h3, FPName, myStyle, 'Format', 'png');
close(h3);

end


% save stats and histograms
if isfield(handles, 'stats') && ~isempty(handles.stats)
stats = handles.stats{channel};


p = str2double(get(handles.edit7, 'String')); % pixel size in nm
bins = str2double(get(handles.edit25, 'String'));

 % cluster area
if base ~= 4
    area = cat(1, stats.Area) * p^2; % cluster area in nm^2
else
    area = stats(:,2); % cluster area in nm^2
end
    h1=figure('visible', 'off');
    hist(area, bins);
    mea = num2str(mean(area));
    med = num2str(median(area));
    st = num2str(std(area));
if channel == 1 %red eventlist
FPName=[dir 'cluster_area_red_' mea '_' med '_' st '.png'];
elseif channel == 2 %green eventlist
FPName=[dir 'cluster_area_green_' mea '_' med '_' st '.png'];
end
    xlabel('cluster area, nm^2')
    ylabel('occurrence')
    hgexport(h1, FPName, hgexport('factorystyle'), 'Format', 'png');
    close(h1);
    
 % diameter
if base ~= 4
    diam = cat(1, stats.EquivDiameter) * p; % cluster diameter in nm
else
    diam = 2 * sqrt(stats(:,2)/pi);
end
    h1=figure('visible', 'off');
    hist(diam, bins);
    mea = num2str(mean(diam));
    med = num2str(median(diam));
    st = num2str(std(diam));
if channel == 1 %red eventlist
FPName=[dir 'cluster_diameter_red_' mea '_' med '_' st '.png'];
elseif channel == 2 %green eventlist
FPName=[dir 'cluster_diameter_green_' mea '_' med '_' st '.png'];
end
    xlabel('cluster diameter, nm')
    ylabel('occurrence')
    hgexport(h1, FPName, hgexport('factorystyle'), 'Format', 'png');
    close(h1);
    
 % eccentricity
 if base ~= 4
    eccentr = cat(1, stats.Eccentricity); % cluster eccentricity 
 else
     eccentr = stats(:,3);
 end
    h1=figure('visible', 'off');
    hist(eccentr, bins);
    mea = num2str(mean(eccentr));
    med = num2str(median(eccentr));
    st = num2str(std(eccentr));
if channel == 1 %red eventlist
FPName=[dir 'cluster_eccentricity_red_' mea '_' med '_' st '.png'];
elseif channel == 2 %green eventlist
FPName=[dir 'cluster_eccentricity_green_' mea '_' med '_' st '.png'];
end
    xlabel('cluster eccentricity')
    ylabel('occurrence')
    hgexport(h1, FPName, hgexport('factorystyle'), 'Format', 'png');
    close(h1);
    
 % events
 if base ~= 4
    I0 = handles.I0{channel};
    s = size(stats, 1);
    NBevents = zeros(s,1);
    for i = 1:s
        ind = stats(i).PixelIdxList;
        NBevents(i) = sum(I0(ind));
    end
 else
     NBevents = stats(:,1);
 end
    h1=figure('visible', 'off');
    hist(NBevents, bins);
    mea = num2str(mean(NBevents));
    med = num2str(median(NBevents));
    st = num2str(std(NBevents));
if channel == 1 %red eventlist
FPName=[dir 'events_in_cluster_red_' mea '_' med '_' st '.png'];
elseif channel == 2 %green eventlist
FPName=[dir 'events_in_cluster_green_' mea '_' med '_' st '.png'];
end
    xlabel('number of events in cluster')
    ylabel('occurrence')
    hgexport(h1, FPName, hgexport('factorystyle'), 'Format', 'png');
    close(h1);
    
 % density
 if base ~= 4
    Dens = zeros(s,1);
    for i = 1:s
        ind = stats(i).PixelIdxList;
        NBeventsi = sum(I0(ind));
        Dens(i) = NBeventsi/(area(i) * 0.001^2); % events/µm^2
    end 
 else
     Dens = 1000 * stats(:,1) ./ stats(:,2); % events/µm^2
 end
    h1=figure('visible', 'off');
    hist(Dens, bins);
    mea = num2str(mean(Dens));
    med = num2str(median(Dens));
    st = num2str(std(Dens));
if channel == 1 %red eventlist
FPName=[dir 'density_in_cluster_red_' mea '_' med '_' st '.png'];
elseif channel == 2 %green eventlist
FPName=[dir 'density_in_cluster_green_' mea '_' med '_' st '.png'];
end
    xlabel('density of events in cluster, 1/µm')
    ylabel('occurrence')
    hgexport(h1, FPName, hgexport('factorystyle'), 'Format', 'png');
    close(h1);
    
clusters = [area, diam, eccentr, NBevents, Dens];
if channel == 1 %red eventlist
FPName=[dir 'clusterStat_red.ascii'];
elseif channel == 2 %green eventlist
FPName=[dir 'clusterStat_green.ascii'];
end    
dlmwrite(FPName, clusters);
end

catch errorObj
    errordlg(getReport(errorObj,'extended','hyperlinks','off'),'Error');
end
    
    

% --------------------------------------------------------------------
function LAS_AF_Callback(hObject, eventdata, handles)
% hObject    handle to LAS_AF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
channel = get(handles.popupmenu1, 'Value');
text{1} = '642-nm eventlist';
text{2} = '488/532-nm eventlist';
if isfield(handles, 'folder')
    pathname = handles.folder;
    [filename, pathname] = uigetfile({'*.ascii'; '*.txt'; '*.ascii'}, text{channel}, pathname);
else
[filename, pathname] = uigetfile({'*.ascii'; '*.txt'; '*.ascii'}, text{channel});
end
if pathname ~= 0
A = importdata([pathname filename]);

if isstruct(A)
A = A.data;
end

% adjust the format to the standard
A = adjustformat(A, 1);

handles.folder = pathname;
handles.AB{channel} = A;
handles.ABmasked{channel} = A;
end
handles.BW = cell(2,1);

%clear old data
cleardata(handles);
guidata(hObject,handles);
updateSRIm(hObject,handles);


% --------------------------------------------------------------------
function QuickPALM_Callback(hObject, eventdata, handles)
% hObject    handle to QuickPALM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
channel = get(handles.popupmenu1, 'Value');
text{1} = '642-nm eventlist';
text{2} = '488/532-nm eventlist';
if isfield(handles, 'folder')
    pathname = handles.folder;
    [filename, pathname] = uigetfile({'*.ascii'; '*.txt'; '*.ascii'}, text{channel}, pathname);
else
[filename, pathname] = uigetfile({'*.ascii'; '*.txt'; '*.ascii'}, text{channel});
end
if pathname ~= 0
A = importdata([pathname filename]);

if isstruct(A)
A = A.data;
end

% adjust the format to the standard
A = adjustformat(A, 2);

handles.folder = pathname;
handles.AB{channel} = A;
handles.ABmasked{channel} = A;
end
handles.BW = cell(2,1);

%clear old data
cleardata(handles);
guidata(hObject,handles);
updateSRIm(hObject,handles);


% --------------------------------------------------------------------
function RapidSTORM_Callback(hObject, eventdata, handles)
% hObject    handle to RapidSTORM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
channel = get(handles.popupmenu1, 'Value');
text{1} = '642-nm eventlist';
text{2} = '488/532-nm eventlist';
if isfield(handles, 'folder')
    pathname = handles.folder;
    [filename, pathname] = uigetfile({'*.ascii'; '*.txt'; '*.ascii'}, text{channel}, pathname);
else
[filename, pathname] = uigetfile({'*.ascii'; '*.txt'; '*.ascii'}, text{channel});
end
if pathname ~= 0
A = importdata([pathname filename]);

if isstruct(A)
A = A.data;
end

% adjust the format to the standard
A = adjustformat(A, 3);

handles.folder = pathname;
handles.AB{channel} = A;
handles.ABmasked{channel} = A;
end
handles.BW = cell(2,1);

%clear old data
cleardata(handles);
guidata(hObject,handles);
updateSRIm(hObject,handles);


% --------------------------------------------------------------------
function Save_mask_Callback(hObject, eventdata, handles)
% hObject    handle to Save_mask (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
channel = get(handles.popupmenu1, 'Value');
BW = handles.BW{channel};
if ~isempty(BW)
pathname = handles.folder;
[FileName,PathName] = uiputfile({'*.tif'; '*.png'; '*.bmp'}, 'Binary mask', pathname);
if FileName ~= 0
FPName = [PathName FileName];
imwrite(BW, FPName);
end
end

function cleardata(handles)
%clear old data

channel = get(handles.popupmenu1, 'Value');
if isfield(handles, 'prevRGB') % image at axes1 (preview)
    handles.prevRGB{channel} = [];
end
if isfield(handles, 'Ripley')
    handles.Ripley(:, channel) = [];
end
if isfield(handles, 'stats')
    handles.stats{channel} = [];
end
if isfield(handles, 'L')
    handles.L(channel, :) = [];
end
if isfield(handles, 'Lim')
    handles.Lim(channel, :) = [];
end
if isfield(handles, 'LimRGB')
    handles.LimRGB(channel, :) = [];
end
cla(handles.axes1); 
cla(handles.axes2); 
cla(handles.axes3); 
cla(handles.axes4); 
cla(handles.axes5); 


% --- Executes when selected object is changed in uipanel8.
function uipanel8_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel8 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
updateLimRGB(hObject, handles);


% --- Executes during object creation, after setting all properties.
function axes1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes1


% --- Executes during object creation, after setting all properties.
function figure1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function uipanel8_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uipanel8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --------------------------------------------------------------------
function uipanel8_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to uipanel8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in popupmenu10.
function popupmenu10_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
setVoronoiThreshold(hObject, handles);
if get(handles.popupmenu10, 'Value') == 1 % L(R)
    set(handles.edit17, 'Enable', 'off'); %brightness
    set(handles.edit7, 'Enable', 'on'); %pixel size
    set(handles.edit8, 'Enable', 'on'); %radius
    set(handles.edit12, 'Enable', 'on'); %noise
elseif get(handles.popupmenu10, 'Value') == 2 % g(r)
    set(handles.edit17, 'Enable', 'off'); %brightness
    set(handles.edit7, 'Enable', 'on'); %pixel size
    set(handles.edit8, 'Enable', 'on'); %radius
    set(handles.edit12, 'Enable', 'on'); %noise
elseif get(handles.popupmenu10, 'Value') == 3 % Gaussian
    set(handles.edit17, 'Enable', 'off'); %brightness
    set(handles.edit7, 'Enable', 'on'); %pixel size
    set(handles.edit8, 'Enable', 'on'); %radius
    set(handles.edit12, 'Enable', 'off'); %noise
elseif get(handles.popupmenu10, 'Value') == 4 % VorDiagram
    set(handles.edit17, 'Enable', 'off'); %brightness
    set(handles.edit7, 'Enable', 'off'); %pixel size
    set(handles.edit8, 'Enable', 'off'); %radius
    set(handles.edit12, 'Enable', 'off'); %noise
elseif get(handles.popupmenu10, 'Value') == 5 % VorDensity
    set(handles.edit17, 'Enable', 'on'); %brightness
    set(handles.edit7, 'Enable', 'on'); %pixel size
    set(handles.edit8, 'Enable', 'off'); %radius
    set(handles.edit12, 'Enable', 'off'); %noise
end



function setVoronoiThreshold(~, handles)
if get(handles.popupmenu10, 'Value') == 4 || get(handles.popupmenu10, 'Value') == 5
    channel = get(handles.popupmenu1, 'Value');
    VoronoiStats = handles.VoronoiStats;
    intersection = VoronoiStats(2,:);
    if ~isempty(intersection{channel})
        set(handles.edit9, 'String', num2str(round(intersection{channel}(1))));
    end
end
if get(handles.popupmenu10, 'Value') == 4
    set(handles.checkbox1, 'Enable', 'off');
    set(handles.edit13, 'Enable', 'off');
else
    set(handles.checkbox1, 'Enable', 'on');
    set(handles.edit13, 'Enable', 'on');
end


% --- Executes during object creation, after setting all properties.
function popupmenu10_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over popupmenu10.
function popupmenu10_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
updateLimRGB(hObject, handles);



function edit17_Callback(hObject, eventdata, handles)
% hObject    handle to edit17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit17 as text
%        str2double(get(hObject,'String')) returns contents of edit17 as a double


% --- Executes during object creation, after setting all properties.
function edit17_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox2.
function checkbox2_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox2


% --- Executes on button press in checkbox3.
function checkbox3_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox3



function edit18_Callback(hObject, eventdata, handles)
% hObject    handle to edit18 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit18 as text
%        str2double(get(hObject,'String')) returns contents of edit18 as a double


% --- Executes during object creation, after setting all properties.
function edit18_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit18 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function microManager_Callback(hObject, eventdata, handles)
% hObject    handle to microManager (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
channel = get(handles.popupmenu1, 'Value');
text{1} = '642-nm eventlist';
text{2} = '488/532-nm eventlist';
if isfield(handles, 'folder')
    pathname = handles.folder;
    [filename, pathname] = uigetfile({'*.ascii'; '*.txt'; '*.ascii'}, text{channel}, pathname);
else
[filename, pathname] = uigetfile({'*.ascii'; '*.txt'; '*.ascii'}, text{channel});
end
if pathname ~= 0
A = importdata([pathname filename]);

if isstruct(A)
A = A.data;
end

% adjust the format to the standard
A = adjustformat(A, 4);

handles.folder = pathname;
handles.AB{channel} = A;
handles.ABmasked{channel} = A;
end
handles.BW = cell(2,1);

%clear old data
cleardata(handles);
guidata(hObject,handles);
updateSRIm(hObject,handles);

% --------------------------------------------------------------------
function SharpViSu_Callback(hObject, eventdata, handles)
% hObject    handle to SharpViSu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
channel = get(handles.popupmenu1, 'Value');
text{1} = 'red eventlist';
text{2} = 'green eventlist';
if isfield(handles, 'folder')
    pathname = handles.folder;
    [filename, pathname] = uigetfile({'*.ascii'; '*.txt'; '*.ascii'}, text{channel}, pathname);
else
[filename, pathname] = uigetfile({'*.ascii'; '*.txt'; '*.ascii'}, text{channel});
end
if pathname ~= 0
A = importdata([pathname filename]);

if isstruct(A)
A = A.data;
end

handles.folder = pathname;
handles.AB{channel} = A;
handles.ABmasked{channel} = A;
end
handles.BW = cell(2,1);

%clear old data
cleardata(handles);
guidata(hObject,handles);
updateSRIm(hObject,handles);
