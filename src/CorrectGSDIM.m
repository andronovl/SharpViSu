function varargout = CorrectGSDIM(varargin)
% CorrectGSDIM MATLAB code for CorrectGSDIM.fig
%      CorrectGSDIM, by itself, creates a new CorrectGSDIM or raises the existing
%      singleton*.
%
%      H = CorrectGSDIM returns the handle to a new CorrectGSDIM or the handle to
%      the existing singleton*.
%
%      CorrectGSDIM('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CorrectGSDIM.M with the given input arguments.
%
%      CorrectGSDIM('Property','Value',...) creates a new CorrectGSDIM or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before CorrectGSDIM_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to CorrectGSDIM_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help CorrectGSDIM

% Last Modified by GUIDE v2.5 21-Jan-2015 15:32:14

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @CorrectGSDIM_OpeningFcn, ...
                   'gui_OutputFcn',  @CorrectGSDIM_OutputFcn, ...
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


% --- Executes just before CorrectGSDIM is made visible.
function CorrectGSDIM_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to CorrectGSDIM (see VARARGIN)

% Choose default command line output for CorrectGSDIM
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes CorrectGSDIM wait for user response (see UIRESUME)
% uiwait(handles.figure1);

%if isempty(gcp('nocreate'))
%    parpool local;
%end
handles.AB = cell(1,2);
guidata(hObject,handles);


% --- Outputs from this function are returned to the command line.
function varargout = CorrectGSDIM_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton13.
function A = pushbutton13_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
format = get(handles.popupmenu21, 'Value');
    extension = '*.*';
if format == 1 % set eventlist extension depending on format
    extension = '*.ascii';
end
if isfield(handles, 'folder')
    pathname = handles.folder;
    [filenameA, pathnameA] = uigetfile({extension}, '642 eventlist', pathname);
else
[filenameA, pathnameA] = uigetfile({extension}, '642 eventlist');
end
if pathnameA ~= 0
A = importdata([pathnameA filenameA]);
pfname = [pathnameA filenameA];
%cd(pathnameB);
set (handles.edit4, 'String' , pfname);

if isstruct(A)
A = A.data;
end

if format ~= 4
A = adjustformat(A, format);    % adjust the format to the standard
end    

handles.folder = pathnameA;
handles.AB{1} = A;
else
    pfname = [];
    set (handles.edit4, 'String' , pfname);
    handles.AB{1} = [];
end

%clear old data
if isfield(handles, 'RGB0')
handles = rmfield(handles,'RGB0');
end
if isfield(handles, 'RGB')
handles = rmfield(handles,'RGB');
end
if isfield(handles, 'ABcorr')
handles = rmfield(handles,'ABcorr');
end
if isfield(handles, 'Res')
handles = rmfield(handles,'Res');
end
if isfield(handles, 'Resd')
handles = rmfield(handles,'Resd');
end
if isfield(handles, 'CA')
handles = rmfield(handles,'CA');
end
if isfield(handles, 'neighbors')
handles = rmfield(handles,'neighbors');
end
if isfield(handles, 'Icoloc')
handles = rmfield(handles,'Icoloc');
end
if isfield(handles, 'offsetAB')
handles = rmfield(handles,'offsetAB');
end
if isfield(handles, 'Image3D')
handles = rmfield(handles,'Image3D');
end
if isfield(handles, 'FRC')
handles = rmfield(handles,'FRC');
end
if isfield(handles, 'driftimageA')
handles = rmfield(handles,'driftimageA');
end
if isfield(handles, 'driftimageB')
handles = rmfield(handles,'driftimageB');
end
if isfield(handles, 'BW')
handles = rmfield(handles,'BW');
end
if isfield(handles, 'Ripley')
handles = rmfield(handles,'Ripley');
end
cla(handles.axes2); %clear the drift axes
cla(handles.axes3); %clear the drift axes
handles.FOV = FOV(handles.AB);
set(handles.edit56, 'String', handles.FOV/1000);
guidata(hObject,handles);



% --- Executes on button press in pushbutton14.
function pushbutton14_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
format = get(handles.popupmenu21, 'Value');
    extension = '*.*';
if format == 1 % set eventlist extension depending on format
    extension = '*.ascii';
end
if isfield(handles, 'folder')
    pathname = handles.folder;
    [filenameB, pathnameB] = uigetfile({extension}, '488 eventlist', pathname);
else
    [filenameB, pathnameB] = uigetfile({extension}, '488 eventlist');
end
if filenameB ~= 0
B = importdata([pathnameB filenameB]); % import eventlist 488
pfname = [pathnameB filenameB];
set (handles.edit3, 'String' , pfname);
handles.folder = pathnameB;

if isstruct(B)
B = B.data;
end

B = adjustformat(B, format);% adjust the format to the standard

handles.AB{2} = B;
else
    pfname = [];
    set (handles.edit3, 'String' , pfname);
    handles.AB{2} = [];
end

%clear old data
if isfield(handles, 'RGB0')
handles = rmfield(handles,'RGB0');
end
if isfield(handles, 'RGB')
handles = rmfield(handles,'RGB');
end
if isfield(handles, 'ABcorr')
handles = rmfield(handles,'ABcorr');
end
if isfield(handles, 'Res')
handles = rmfield(handles,'Res');
end
if isfield(handles, 'Resd')
handles = rmfield(handles,'Resd');
end
if isfield(handles, 'CA')
handles = rmfield(handles,'CA');
end
if isfield(handles, 'neighbors')
handles = rmfield(handles,'neighbors');
end
if isfield(handles, 'Icoloc')
handles = rmfield(handles,'Icoloc');
end
if isfield(handles, 'offsetAB')
handles = rmfield(handles,'offsetAB');
end
if isfield(handles, 'Image3D')
handles = rmfield(handles,'Image3D');
end
if isfield(handles, 'FRC')
handles = rmfield(handles,'FRC');
end
if isfield(handles, 'driftimageA')
handles = rmfield(handles,'driftimageA');
end
if isfield(handles, 'driftimageB')
handles = rmfield(handles,'driftimageB');
end
if isfield(handles, 'BW')
handles = rmfield(handles,'BW');
end
if isfield(handles, 'Ripley')
handles = rmfield(handles,'Ripley');
end
cla(handles.axes2); %clear the drift axes
cla(handles.axes3); %clear the drift axes
handles.FOV = FOV(handles.AB);
set(handles.edit56, 'String', handles.FOV/1000);
guidata(hObject,handles);


% --- Executes on button press in pushbutton17.
function [handles] = pushbutton17_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

stringsteps642 = get(handles.edit28, 'string');
s(1) = str2double(stringsteps642);
stringsteps488 = get(handles.edit29, 'string');
s(2) = str2double(stringsteps488);
corrpixsize(1) = str2double(get(handles.edit35, 'string'));
corrpixsize(2) = str2double(get(handles.edit36, 'string'));

offsetAB = cell(2,1);
totaldrift = cell(2,1);
corrected = get(handles.checkbox24, 'Value');%corrected lists
if ~corrected
for i = 1:2
if isfield(handles, 'AB') && ~isempty(handles.AB{i})
    AB = handles.AB{i};
    [offsetAB{i}, totaldrift{i}] = offset(AB, corrpixsize(i), s(i));
else
    offsetAB{i} = [NaN NaN];
    totaldrift{i} = [NaN NaN];
end
end
elseif corrected
for i = 1:2
if isfield(handles, 'ABcorr') && ~isempty(handles.ABcorr{i})
    AB = handles.ABcorr{i};
    [offsetAB{i}, totaldrift{i}] = offset(AB, corrpixsize(i), s(i));
else
    offsetAB{i} = [NaN NaN];
    totaldrift{i} = [NaN NaN];
end
end
end
% [1st_arrow_begin_x 2nd_arrow_begin_x 3d_arrow_begin_x ...], 
% [1st_arrow_begin_y 2nd_arrow_begin_y 3d_arrow_begin_y ...], 
% [1st_arrow_length_x 2nd_arrow_length_x 3d_arrow_length_x ...],
% [1st_arrow_length_y 2nd_arrow_length_y 3d_arrow_length_y ...]
% nth_arrow_length_x(y) = offset1[n, 1(2)] from offset.m
% nth_arrow_begin_x(y) = sum of the lengths of the arrows until the (n-1)th
cla(handles.axes2);
if ~isnan(totaldrift{1})
    begA = [0 0; cumsum(totaldrift{1}(1:end-1, 1)) cumsum(totaldrift{1}(1:end-1, 2))];
    lenxA = totaldrift{1}(:,1);
    lenyA = totaldrift{1}(:,2);
    mod = zeros(size(totaldrift{1},1), 1);
    for i = 1 : size(totaldrift{1},1)
    mod(i) = norm(totaldrift{1}(i,:));
    end
    handles.driftimageA = quiver('Parent', handles.axes2, begA(:,1), begA(:,2), lenxA, lenyA, 0, 'LineWidth', 1.5, 'MaxHeadSize', 0.5, 'Color', [1 0 0]);
    set(handles.axes2,'xaxislocation','top','ydir','reverse');
    
    if get(handles.checkbox56, 'Value')% save the drift picture
    pathname = handles.folder;
    i = 0;
    while exist([pathname 'drift_red_' num2str(i) '.png'], 'file')
        i = i + 1;
    end
    FPName = [pathname 'drift_red_' num2str(sum(mod)) '.png'];
    h2 = figure('visible', 'off');
    quiver(begA(:,1), begA(:,2), lenxA, lenyA, 0, 'LineWidth', 3, 'MaxHeadSize', 0.5, 'Color', [1 0 0]);
    xlabel('drift X, nm', 'FontSize', 18);
    ylabel('drift Y, nm', 'FontSize', 18);
    set(gca, 'xaxislocation', 'top', 'ydir', 'reverse', 'FontSize', 18, 'LineWidth', 1);
    %axis ([-50 10 -35 25]);
    myStyle = hgexport('factorystyle');
    myStyle.Format = 'png';
    myStyle.Width = 6;
    myStyle.Height = 6;
    myStyle.Resolution = 300;
    myStyle.Units = 'inch';
    myStyle.FixedFontSize = 8;
    hgexport(h2, FPName, myStyle, 'Format', 'png');
    close(h2);
    end
end
cla(handles.axes3);
if ~isnan(totaldrift{2})
    begB = [0 0; cumsum(totaldrift{2}(1:end-1, 1)) cumsum(totaldrift{2}(1:end-1, 2))];
    lenxB = totaldrift{2}(:,1);
    lenyB = totaldrift{2}(:,2);
    mod = zeros(size(totaldrift{2},1), 1);
    for i = 1 : size(totaldrift{2},1)
    mod(i) = norm(totaldrift{2}(i,:));
    end
    handles.driftimageB = quiver('Parent', handles.axes3, begB(:,1), begB(:,2), lenxB, lenyB, 0, 'LineWidth', 1.5, 'MaxHeadSize', 0.5, 'Color', [0 0.5 0]);
    set(handles.axes3,'xaxislocation','top','ydir','reverse');
    
    if get(handles.checkbox56, 'Value')
    % save the drift picture
    pathname = handles.folder;
    i = 0;
    while exist([pathname 'drift_green_' num2str(i) '.png'], 'file')
        i = i + 1;
    end
    FPName = [pathname 'drift_green_' num2str(sum(mod)) '.png'];
    h2 = figure('visible', 'off');
    quiver(begB(:,1), begB(:,2), lenxB, lenyB, 0, 'LineWidth', 3, 'MaxHeadSize', 0.5, 'Color', [0 0.5 0]);
    xlabel('drift X, nm', 'FontSize', 18);
    ylabel('drift Y, nm', 'FontSize', 18);
    set(gca, 'xaxislocation', 'top', 'ydir', 'reverse', 'FontSize', 18, 'LineWidth', 1);
    myStyle = hgexport('factorystyle');
    myStyle.Format = 'png';
    myStyle.Width = 6;
    myStyle.Height = 6;
    myStyle.Resolution = 300;
    myStyle.Units = 'inch';
    myStyle.FixedFontSize = 8;
    hgexport(h2, FPName, myStyle, 'Format', 'png');
    close(h2);
    end
end
handles.totaldrift = totaldrift;
handles.offsetAB = offsetAB;
guidata(hObject,handles);




% --- Executes on button press in pushbutton18.
function pushbutton18_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton18 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Ic = handles.Ic;
I0 = handles.I0;
p = handles.FOV/size(Ic{1},1); % pixel size of Ic, I0 in nm

[~, RGB, b] = contrastenh(Ic, p, 'corrected super-resolution image');


if draworig == 1
    [~, RGB0, c] = contrastenh(I0, p, 'non-corrected original super-resolution image');
else
    c = zeros(2,1);
end

evvalue01 = num2str(c(1));
evvalue02 = num2str(c(2));
evvalue1 = num2str(b(1));
evvalue2 = num2str(b(2));
if draworig == 1
saveorig = get(handles.checkbox26, 'Value'); %auto save original image
end
savecorr = get(handles.checkbox27, 'Value'); %auto save corrected image
pathname = handles.folder;
if draworig == 1 && saveorig == 1 
    if histmode == 1
    filenameorig = ['0h_' evvalue01 '_' evvalue02 '.tif'];
    elseif gaussmode == 1
    filenameorig = ['0g_' evvalue01 '_' evvalue02 '.tif'];
    end
    FPNameorig = [pathname filenameorig];
    imwrite(RGB0, FPNameorig);
end
if savecorr == 1
    if histmode == 1
    filenamecorr = ['1h_' evvalue1 '_' evvalue2 '.tif'];
    elseif gaussmode == 1
    filenamecorr = ['1g_' evvalue1 '_' evvalue2 '.tif'];
    end
    FPNamecorr = [pathname filenamecorr];
    imwrite(RGB, FPNamecorr);
end


handles.RGB = RGB;
if (draworig == 1)
handles.RGB0 = RGB0;
end

guidata(hObject,handles);


% --- Executes on button press in pushbutton19.
function pushbutton19_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
RGB = handles.RGB;
RGB0 = handles.RGB0;
pathname = handles.folder;
[FileName,PathName] = uiputfile({'*.tif'}, 'Corrected RGB', pathname);
if FileName ~= 0
FPName=[PathName FileName];
imwrite(RGB, FPName);
end
[FileName0,PathName0] = uiputfile({'*.tif'}, 'Non-corrected RGB', pathname);
if FileName0 ~= 0
FPName0 = [PathName0 FileName0];
imwrite(RGB0, FPName0);
end

% --- Executes on button press in checkbox1.
function checkbox1_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox1


function [offset, offset1] = offset(A, corrpixsize, s)
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
I12{i} = draw(A, corrpixsize, 0, ko, ke);
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
end


 
%drawing a histogram image
function [I] = draw(A, p, fov, ko, ke)
% draws a histogram image I on frames from ko to ke with pixel size p in a
% square with a side = fov
if ~exist('ko', 'var')
    ko = 1;
end
if ~exist('ke', 'var')
    ke = size(A, 1);
end

if ~exist('fov', 'var') || fov == 0
    fov = FOV(A);
end
Anew(:,1) = A(ko:ke,5);
Anew(:,2) = A(ko:ke,4);
edges = 0:p:fov-p;
edges = {edges, edges};
I = hist3 (Anew, 'Edges', edges);


% --- Executes on button press in pushbutton21.
function pushbutton21_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton21 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
contrastenh(handles.AB, handles.ABcorr, handles.folder);



function [Acorr] = corrend(A, offset) 
% shifts the events towards the end of the acquisition using offset as
% the frame-to-frame drift value (for correction of the channel imaged first [642 nm by default])
s = 1 + size(offset,1); % number of parts
l = size(A, 1); % number of rows
border_event = floor(l/s); %the index of the final frame of I1
middle_event = zeros(s,1); %index of the middle event of the eventlist i
middle_frame = zeros(s,1);

for i = 1:s
    ko = 1 + border_event * (i - 1);
    ke = border_event * i;
middle_event(i) = round((ke - ko) / 2 + ko - 1/2);
middle_frame(i) = A(middle_event(i),2);
end

totoffset = zeros(s-1,2);

for i = 1:s-1
offset_length = middle_frame(i+1) - middle_frame(i);
totoffset(i,:) = offset(i,:) * offset_length; %total drift between two consecutive eventlists
end

Acorr = cell(s-1,1);

%cut the entire eventlist into correspondings sublists
if s > 2
    Acorr{s-1} = A(middle_event(s-1)+1:end, :); % the last sublist
else
    Acorr{1} = A; % the only sublist in case if s == 2
end

    ln = size(Acorr{s-1},1);
    last_frame = Acorr{s-1}(ln,2);
    
for n = 1 : ln
for k = 4 : 5
    Acorr{s-1}(n,k) = Acorr{s-1}(n,k) + offset(s-1,k-3) * (last_frame - Acorr{s-1}(n,2)); %correction of the last sublist
end
end

sumoffset = zeros(1,2);
if s > 3
for i = s-2 : -1 : 2
    Acorr{i} = A(middle_event(i)+1:middle_event(i+1), :);
    ln = size(Acorr{i},1);
    last_frame = Acorr{i}(ln,2);
for n = 1 : ln
for k = 4 : 5
    Acorr{i}(n,k) = Acorr{i}(n,k) + 1.5 * totoffset(s-1,k-3) + sumoffset(k-3) + offset(i,k-3) * (last_frame - Acorr{i}(n,2));
end
end
    sumoffset = sumoffset + totoffset(i, :);
end
end

if s > 2
    Acorr{1} = A(1:middle_event(2),:);
    ln = size(Acorr{1},1);
    last_frame = Acorr{1}(ln,2);
for n = 1 : ln
for k = 4 : 5
    Acorr{1}(n,k) = Acorr{1}(n,k) + 1.5 * totoffset(s-1,k-3) + sumoffset(k-3) + offset(1,k-3) * (last_frame - Acorr{1}(n,2));
end
end
end

Acorr = cell2mat(Acorr);

function [Acorr] = corrbeg(A, offset) 
% shifts the events to the beginning of the acquisition using offset as
% the frame-to-frame drift value (for correction of the channel imaged second [488 nm or 532 nm by default])

s = 1 + size(offset,1);
l = size(A, 1); % number of rows
border_event = floor(l/s); %the index of the final frame of I1
middle_event = zeros(s,1); %index of the middle event of the eventlist i
middle_frame = zeros(s,1);
for i = 1:s
    ko = 1 + border_event * (i - 1);
    ke = border_event * i;
middle_event(i) = round((ke - ko) / 2 + ko - 1/2);
middle_frame(i) = A(middle_event(i),2);
end
totoffset = zeros(s-1,2);
for i = 1:s-1
offset_length = middle_frame(i+1) - middle_frame(i);
totoffset(i,:) = offset(i,:) * offset_length; %total drift between two consecutive eventlists
end

Acorr = cell(s-1,1);
if s > 2
    Acorr{1} = A(1:middle_event(2),:);
else
    Acorr{1} = A;
end
    ln = size(Acorr{1},1);
for n = 1 : ln
for k = 4 : 5
    Acorr{1}(n,k) = Acorr{1}(n,k) - offset(1,k-3) * Acorr{1}(n,2);
end
end

sumoffset = zeros(1,2);
if s > 3
for i = 2:s-2
    Acorr{i} = A(middle_event(i)+1:middle_event(i+1), :);
    ln = size(Acorr{i},1);
    first_frame = Acorr{i}(1,2);
for n = 1 : ln
for k = 4 : 5
    Acorr{i}(n,k) = Acorr{i}(n,k) - 1.5 * totoffset(1,k-3) - sumoffset(k-3) - offset(i,k-3) * (Acorr{i}(n,2) - first_frame);
end
end
    sumoffset = sumoffset + totoffset(i, :);
end
end

if s > 2
    Acorr{s-1} = A(middle_event(s-1)+1:end, :);
    ln = size(Acorr{s-1},1);
    first_frame = Acorr{s-1}(1,2);
for n = 1 : ln
for k = 4 : 5
    Acorr{s-1}(n,k) = Acorr{s-1}(n,k) - 1.5 * totoffset(1,k-3) - sumoffset(k-3) - offset(s-1,k-3) * (Acorr{s-1}(n,2) - first_frame);
end
end
end

Acorr = cell2mat(Acorr);


% --- Executes on button press in pushbutton43.
function pushbutton43_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton43 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
A = handles.ABcorr{1};
pathname = handles.folder;
[FileName,PathName] = uiputfile({'*.ascii'}, '647 corrected eventlist', pathname);
if FileName ~= 0
FPName=[PathName FileName];
dlmwrite(FPName, A);
end

% --- Executes on button press in pushbutton42.
function pushbutton42_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton42 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
A = handles.ABcorr{2};
pathname = handles.folder;
[FileName,PathName] = uiputfile({'*.ascii'}, '488 corrected eventlist', pathname);
if FileName ~= 0
FPName=[PathName FileName];
dlmwrite(FPName, A);
end


% --- Executes on button press in pushbutton56.
function pushbutton56_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton56 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
AB = handles.AB;

chroma = get(handles.checkbox21, 'Value'); %correction of chroma aberrations

% correction of chromatic aberrations
    el488 = get(handles.radiobutton11, 'Value');
    el532 = get(handles.radiobutton10, 'Value');
if ~isempty(AB{2}) && chroma && el488
	AB{2} = shift488(AB{2});
elseif ~isempty(AB{2}) && chroma && el532
	AB{2} = shift532(AB{2});
end

iter(1) = str2double(get(handles.edit61, 'String'));%642 drift iterations
iter(2) = str2double(get(handles.edit62, 'String'));%488 drift iterations

%correction of drift
offsetAB = cell(2,1);
ABcorr = AB;
if max(iter) > 0
checkbox24old = get(handles.checkbox24, 'Value');
h = waitbar(0, 'Iterative drift correction...');
for j = 1:max(iter)
if j > 1
    ABcorr = handles.ABcorr;
end
for i = 1:2    
    
if (iter(i) < 1) || (iter(i) - j < 0) || ~isfield(handles, 'offsetAB')
    offsetAB{i} = [NaN NaN];
else
    offsetAB{i} = handles.offsetAB{i};
end

if ~isnan(offsetAB{i})
    if i == 1
    ABcorr{i} = corrend(ABcorr{i}, offsetAB{i});
    assignin('base', ['dr' num2str(j)] ,handles.totaldrift{i});
    elseif i == 2
    ABcorr{i} = corrbeg(ABcorr{i}, offsetAB{i});
    end
end

end
handles.ABcorr = ABcorr;
set(handles.checkbox24, 'Value', 1);
handles = pushbutton17_Callback(hObject, eventdata, handles);
waitbar(j/max(iter));
end
close(h);

set(handles.checkbox24, 'Value', checkbox24old);
end

save = get(handles.checkbox25, 'Value'); %auto save corrected eventlists
if save == 1
    el488 = get(handles.radiobutton11, 'Value');
    el532 = get(handles.radiobutton10, 'Value');
    pathname = handles.folder;
    if el488 == 1
    filename488 = 'el_488_corrected.ascii';
    elseif el532 == 1
    filename488 = 'el_532_corrected.ascii';
    end
    FPName488 = [pathname filename488];
    filename642 = 'el_642_corrected.ascii';
    FPName642 = [pathname filename642];
    dlmwrite(FPName488, ABcorr{2});
    dlmwrite(FPName642, ABcorr{1});
end

handles.ABcorr = ABcorr;
guidata(hObject,handles);


% --- Executes on button press in pushbutton55.
function pushbutton55_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton55 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Image3D = handles.Image3D;
l = size(Image3D, 2);
pathname = handles.folder;
[FileName,PathName] = uiputfile({'*.tif'}, 'Choose a folder for 3D-stack', pathname);
if FileName ~= 0
for i = 1:l
filename = num2str(i);
filename = strcat(filename, '.tif');
FPName=[PathName filename];
imwrite(Image3D{i}, FPName);
end
end

% --- Executes on button press in pushbutton53.
function pushbutton53_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton53 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
oil = get(handles.radiobutton9, 'Value');

if oil == 1
    m = 1;
else
    m = 0.79;
end

%determine the color of the second eventlist
el488 = get(handles.radiobutton11, 'Value');
el532 = get(handles.radiobutton10, 'Value');

AB = handles.ABcorr;
        
if ~isempty (AB{1})
    AB{1} = zcalc(AB{1}, 1, m);
end
if ~isempty (AB{2})
    if el488
    AB{2} = zcalc(AB{2}, 3, m);
    elseif el532
    AB{2} = zcalc(AB{2}, 2, m);
    end
end

handles.ABcorr = AB;
guidata(hObject,handles);



% --- Executes on button press in pushbutton57.
function pushbutton57_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton57 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ABcorr = handles.ABcorr;
pixsize = str2double(get(handles.edit42, 'string'));

fov = FOV(ABcorr);
I = cell(1,2);
for i = 1:2
if ~isempty(ABcorr{i})
    I{i} = draw(ABcorr{i}, pixsize, fov);
else
    I{i} = zeros(fov/pixsize);
end
end

RGB = cat(3, I{1}, I{2}, zeros(fov/pixsize));
figure;
BW = roipoly(RGB);
if ~isempty(BW)
close;
if ~isempty(ABcorr{1})
ABcorr{1} = parroifilter( ABcorr{1}, BW, 1);
end
if ~isempty(ABcorr{2})
ABcorr{2} = parroifilter( ABcorr{2}, BW, 1);
end
save = get(handles.checkbox25, 'Value'); %auto save corrected eventlists
if save == 1
    el488 = get(handles.radiobutton11, 'Value');
    el532 = get(handles.radiobutton10, 'Value');
    pathname = handles.folder;
    if el488 == 1
    filename488 = 'el_488_corrected.ascii';
    elseif el532 == 1
    filename488 = 'el_532_corrected.ascii';
    end
    FPName488 = [pathname filename488];
    filename642 = 'el_642_corrected.ascii';
    FPName642 = [pathname filename642];
    dlmwrite(FPName488, ABcorr{2});
    dlmwrite(FPName642, ABcorr{1});
end
handles.ABcorr = ABcorr;
handles.BW = BW;
guidata(hObject,handles);
end


% --- Executes on button press in pushbutton58.
function pushbutton58_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton58 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ABcorr = handles.ABcorr;
pixsize = str2double(get(handles.edit42, 'string'));


fov = FOV(ABcorr);
I = cell(2,1);
for i = 1:2
if (isempty(ABcorr{i}) == 0)
    I{i} = draw(ABcorr{i}, pixsize, fov);
else
    I{i} = zeros(fov/pixsize);
end
end

RGB = cat(3, I{1}, I{2}, zeros(fov/pixsize));
figure;
BW = roipoly(RGB);
if ~isempty(BW)
close;
for i = 1:2
    if ~isempty(ABcorr{i})
ABcorr{i} = parroifilter( ABcorr{i}, BW, 0);
    end
end
save = get(handles.checkbox25, 'Value'); %auto save corrected eventlists
if save == 1
    el488 = get(handles.radiobutton11, 'Value');
    el532 = get(handles.radiobutton10, 'Value');
    pathname = handles.folder;
    if el488 == 1
    filename488 = 'el_488_corrected.ascii';
    elseif el532 == 1
    filename488 = 'el_532_corrected.ascii';
    end
    FPName488 = [pathname filename488];
    filename642 = 'el_642_corrected.ascii';
    FPName642 = [pathname filename642];
    dlmwrite(FPName488, ABcorr{2});
    dlmwrite(FPName642, ABcorr{1});
end
handles.ABcorr = ABcorr;
guidata(hObject,handles);
end


% --- Executes on button press in pushbutton59.
function pushbutton59_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton59 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
stepsstr = get(handles.edit27, 'string');
steps = str2double(stepsstr);
AB = handles.ABcorr;

FRC = cell(1,2);
resol = cell(1,2);
%parfor i = 1:2    
for i = 1:2  
if ~isempty (AB{i})
    [FRC{i}, resol{i}] = resolution(AB{i}, steps);
end
end
set (handles.text83, 'String' , resol{2});
set (handles.text81, 'String' , resol{1});


save = get(handles.checkbox37, 'Value'); %auto save FRC curve
if save == 1
    el488 = get(handles.radiobutton11, 'Value');
    el532 = get(handles.radiobutton10, 'Value');
    pathname = handles.folder;
    for i = 1:2
        if i == 1
           FPName = [pathname 'FRC647.png'];
        elseif el488 == 1
           FPName = [pathname 'FRC488.png'];
        elseif el532 == 1
           FPName = [pathname 'FRC532.png'];
        end
    h2 = figure('visible', 'off');
    if ~isempty (FRC{i})
    FRC1 = FRC{i};
    plot(FRC1(:,1),FRC1(:,3),'-',FRC1(:,1),FRC1(:,2),'.', 'LineWidth', 1, 'MarkerSize', 8);
    xlabel('Spatial frequency (nm^{-1})');
    ylabel('FRC');
    set(gca, 'XTick', 0:0.01:0.05, 'YTick', 0:0.2:1, 'FontSize', 14, 'LineWidth', 1);
    axis([0 0.05 0 1]);
    myStyle = hgexport('factorystyle');
    myStyle.Format = 'png';
    myStyle.Width = 6;
    myStyle.Height = 4;
    myStyle.Resolution = 300;
    myStyle.Units = 'inch';
    myStyle.FixedFontSize = 8;
    hgexport(h2, FPName, myStyle, 'Format', 'png');
    close(h2);
    end
    end
end

handles.FRC = FRC;
guidata(hObject,handles);

% --- Executes on button press in pushbutton60.
function pushbutton60_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton60 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
FRC = handles.FRC;
if ~isempty (FRC{2})
    FRC1 = FRC{2};
    figure(2);
    plot(FRC1(:,1),FRC1(:,3),'-',FRC1(:,1),FRC1(:,2),'.', 'LineWidth', 1, 'MarkerSize', 8);
    xlabel('Spatial frequency, (nm^{-1})');
    ylabel('FRC');
    set(gca, 'XTick', 0:0.01:0.05, 'YTick', 0:0.2:1, 'FontSize', 14, 'LineWidth', 1);
    axis([0 0.05 0 1]);
end


% --- Executes on button press in pushbutton62.
function pushbutton62_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton62 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
FRC = handles.FRC;
if ~isempty (FRC{1})
    FRC1 = FRC{1};
    figure(1);
    plot(FRC1(:,1),FRC1(:,3),'-',FRC1(:,1),FRC1(:,2),'.', 'LineWidth', 1, 'MarkerSize', 8);
    xlabel('Spatial frequency, (nm^{-1})');
    ylabel('FRC');
    set(gca, 'XTick', 0:0.01:0.05, 'YTick', 0:0.2:1, 'FontSize', 14, 'LineWidth', 1);
    axis([0 0.05 0 1]);
end


% --- Executes on button press in pushbutton63.
function pushbutton63_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton63 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
timeincolor(handles.ABcorr, handles.AB, handles.folder);



% --- Executes on button press in pushbutton77.
function pushbutton77_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton77 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
AB = handles.ABcorr;
[~, Image3D] = viewZ(AB);

savestack = get(handles.checkbox28, 'Value'); %auto save z-stack
if savestack == 1
    l = size(Image3D, 2);
    pathname = handles.folder;
    folder = [pathname 'stack\'];
    mkdir(folder);
    for i = 1:l
        filename = num2str(i);
        filename = strcat(filename, '.tif');
        FPName=[folder filename];
        imwrite(Image3D{i}, FPName, 'Compression', 'lzw');
    end
end
handles.Image3D = Image3D;
guidata(hObject,handles);



% --- Executes on button press in pushbutton78.
function pushbutton78_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton78 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ABcorr = handles.ABcorr;
[~, ABcorr] = filtration(ABcorr);

save = get(handles.checkbox25, 'Value'); %auto save corrected eventlists
if save == 1
    el488 = get(handles.radiobutton11, 'Value');
    el532 = get(handles.radiobutton10, 'Value');
    pathname = handles.folder;
    if el488 == 1
    filename488 = 'el_488_corrected.ascii';
    elseif el532 == 1
    filename488 = 'el_532_corrected.ascii';
    end
    FPName488 = [pathname filename488];
    filename642 = 'el_642_corrected.ascii';
    FPName642 = [pathname filename642];
    dlmwrite(FPName488, ABcorr{2});
    dlmwrite(FPName642, ABcorr{1});
end

handles.ABcorr = ABcorr;
guidata(hObject,handles);


% --- Executes on button press in pushbutton80.
function pushbutton80_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton80 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA).
if isfield (handles, 'folder')
    Zcalibration(handles.folder);
else
    Zcalibration;
end

% --- Executes on button press in pushbutton83.
function pushbutton83_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton83 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Chroma_calibration;


% --- Executes on button press in checkbox56.
function checkbox56_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox56 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox56


% --- Executes on button press in checkbox25.
function checkbox25_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox25 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox25


% --- Executes on button press in checkbox28.
function checkbox28_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox28 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox28


% --- Executes on button press in checkbox37.
function checkbox37_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox37 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox37


% --- Executes on button press in checkbox53.
function checkbox53_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox53 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox53


% --- Executes on button press in checkbox54.
function checkbox54_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox54 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox54



function edit27_Callback(hObject, eventdata, handles)
% hObject    handle to edit27 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit27 as text
%        str2double(get(hObject,'String')) returns contents of edit27 as a double


% --- Executes on button press in checkbox24.
function checkbox24_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox24 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox24


% --- Executes on button press in checkbox21.
function checkbox21_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox21 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox21
