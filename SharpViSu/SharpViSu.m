function varargout = SharpViSu(varargin)
% SharpViSu MATLAB code for SharpViSu.fig
%      SharpViSu, by itself, creates a new SharpViSu or raises the existing
%      singleton*.
%
%      H = SharpViSu returns the handle to a new SharpViSu or the handle to
%      the existing singleton*.
%
%      SharpViSu('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SharpViSu.M with the given input arguments.
%
%      SharpViSu('Property','Value',...) creates a new SharpViSu or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before SharpViSu_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to SharpViSu_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help SharpViSu

% Last Modified by GUIDE v2.5 20-Jan-2016 16:44:56

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @SharpViSu_OpeningFcn, ...
                   'gui_OutputFcn',  @SharpViSu_OutputFcn, ...
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


% --- Executes just before SharpViSu is made visible.
function SharpViSu_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to SharpViSu (see VARARGIN)

% Choose default command line output for SharpViSu
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes SharpViSu wait for user response (see UIRESUME)
% uiwait(handles.figure1);

%if isempty(gcp('nocreate'))
%    parpool local;
%end
handles.AB = cell(1,2);
guidata(hObject,handles);


% --- Outputs from this function are returned to the command line.
function varargout = SharpViSu_OutputFcn(hObject, eventdata, handles) 
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
try

format = get(handles.popupmenu21, 'Value');
    extension = '*.*';
if format == 1 % set eventlist extension depending on format
    extension = '*.ascii';
elseif format == 4 % set eventlist extension depending on format
    extension = '*.txt';
end
if isfield(handles, 'folder')
    pathname = handles.folder;
    [filenameA, pathnameA] = uigetfile({extension}, '642 eventlist', pathname);
else
[filenameA, pathnameA] = uigetfile({extension}, '642 eventlist');
end
if pathnameA ~= 0
    if format == 4 % µManager: convert commas to dots
        file    = memmapfile([pathnameA filenameA], 'writable', true );
        comma   = uint8(',');
        point   = uint8('.');
        file.Data( transpose( file.Data==comma) ) = point;
    end
A = importdata([pathnameA filenameA]);
pfname = [pathnameA filenameA];
%cd(pathnameB);
set (handles.edit4, 'String' , pfname);

if isstruct(A)
A = A.data;
end

if format ~= 5
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
set (handles.text83, 'String' , '—'); %clear the resolution
set (handles.text81, 'String' , '—'); %clear the resolution
handles.FOV = FOV(handles.AB);
set(handles.edit56, 'String', handles.FOV/1000);
guidata(hObject,handles);

catch errorObj
    errordlg(getReport(errorObj,'extended','hyperlinks','off'),'Error');
end



% --- Executes on button press in pushbutton14.
function pushbutton14_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try

format = get(handles.popupmenu21, 'Value');
    extension = '*.*';
if format == 1 % set eventlist extension depending on format
    extension = '*.ascii';
elseif format == 4 % set eventlist extension depending on format
    extension = '*.txt';
end
if isfield(handles, 'folder')
    pathname = handles.folder;
    [filenameB, pathnameB] = uigetfile({extension}, '488 eventlist', pathname);
else
    [filenameB, pathnameB] = uigetfile({extension}, '488 eventlist');
end

if pathnameB ~= 0
    if format == 4 % µManager: convert commas to dots
        file    = memmapfile([pathnameB filenameB], 'writable', true );
        comma   = uint8(',');
        point   = uint8('.');
        file.Data( transpose( file.Data==comma) ) = point;
    end
B = importdata([pathnameB filenameB]);
pfname = [pathnameB filenameB];
%cd(pathnameB);
set (handles.edit3, 'String' , pfname);

if isstruct(B)
B = B.data;
end

if format ~= 5
B = adjustformat(B, format);    % adjust the format to the standard
end    

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
set (handles.text83, 'String' , '—'); %clear the resolution
set (handles.text81, 'String' , '—'); %clear the resolution
handles.FOV = FOV(handles.AB);
set(handles.edit56, 'String', handles.FOV/1000);
guidata(hObject,handles);

catch errorObj
errordlg(getReport(errorObj,'extended','hyperlinks','off'),'Error');
end


% --- Executes on button press in pushbutton17.
function [handles] = pushbutton17_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try

method = get(handles.popupmenu23, 'Value');
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
    [offsetAB{i}, totaldrift{i}] = offset(AB, corrpixsize(i), s(i), method);
else
    offsetAB{i} = [NaN NaN];
    totaldrift{i} = [NaN NaN];
end
end
elseif corrected
for i = 1:2
if isfield(handles, 'ABcorr') && ~isempty(handles.ABcorr{i})
    AB = handles.ABcorr{i};
    [offsetAB{i}, totaldrift{i}] = offset(AB, corrpixsize(i), s(i), method);
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


if isempty(handles.AB{1}) && isempty(handles.AB{2})
    h = errordlg('Please load data first','No data');
end

catch errorObj
errordlg(getReport(errorObj,'extended','hyperlinks','off'),'Error');
end



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



% --- Executes on button press in pushbutton21.
function pushbutton21_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton21 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try
contrastenh(handles.AB, handles.ABcorr, handles.folder);
catch errorObj
errordlg(getReport(errorObj,'extended','hyperlinks','off'),'Error');
end




% --- Executes on button press in pushbutton43.
function pushbutton43_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton43 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try

A = handles.ABcorr{1};
pathname = handles.folder;
[FileName,PathName,FilterIndex] = uiputfile( ...
{'*.ascii', 'SharpViSu (*.ascii)';...
 '*.2d', 'ViSP 2D File (*.2d)';...
 '*.2dlp', 'ViSP 2D Localization Precision File (*.2dlp)';...
 '*.3d', 'ViSP 3D File (*.3d)';...
 '*.3dlp', 'ViSP 3D Localization Precision File (*.3dlp)'},...
 '647 corrected eventlist', pathname);
if FilterIndex ~= 0
    FPName=[PathName FileName];
    if FilterIndex == 1
        dlmwrite(FPName, A, '\t');
    else
        SaveViSP( A, FilterIndex - 1, FPName );
    end
end

catch errorObj
errordlg(getReport(errorObj,'extended','hyperlinks','off'),'Error');
end

% --- Executes on button press in pushbutton42.
function pushbutton42_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton42 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try

A = handles.ABcorr{2};
pathname = handles.folder;
[FileName,PathName,FilterIndex] = uiputfile( ...
{'*.ascii', 'SharpViSu (*.ascii)';...
 '*.2d', 'ViSP 2D File (*.2d)';...
 '*.2dlp', 'ViSP 2D Localization Precision File (*.2dlp)';...
 '*.3d', 'ViSP 3D File (*.3d)';...
 '*.3dlp', 'ViSP 3D Localization Precision File (*.3dlp)'},...
 '488 corrected eventlist', pathname);
if FilterIndex ~= 0
    FPName=[PathName FileName];
    if FilterIndex == 1
        dlmwrite(FPName, A, '\t');
    else
        SaveViSP( A, FilterIndex - 1, FPName );
    end
end

catch errorObj
errordlg(getReport(errorObj,'extended','hyperlinks','off'),'Error');
end


% --- Executes on button press in pushbutton56.
function pushbutton56_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton56 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try

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
cons = get(handles.checkbox58, 'Value'); % consecutively acquired channels

%correction of drift
offsetAB = cell(2,1);
ABcorr = AB;
if max(iter) > 0
checkbox24old = get(handles.checkbox24, 'Value');
h = waitbar(0, 'Iterative drift correction...');

for j = 1 : max(iter)

if j > 1
    ABcorr = handles.ABcorr;
end  

for i = 1:2
if (iter(i) - j < 0) || ~isfield(handles, 'offsetAB')
    offsetAB{i} = [NaN NaN];
else
    offsetAB{i} = handles.offsetAB{i};
end

if ~isnan(offsetAB{i})
    if i == 1
        if cons
            ABcorr{i} = corrend(ABcorr{i}, offsetAB{i});
        else
            ABcorr{i} = corrbeg(ABcorr{i}, offsetAB{i});
        end
    %assignin('base', ['dr' num2str(j)], handles.totaldrift{i});
    elseif i == 2
        ABcorr{i} = corrbeg(ABcorr{i}, offsetAB{i});
    end
end

end
handles.ABcorr = ABcorr;
set(handles.checkbox24, 'Value', 1);
if ~isnan(offsetAB{1}(1,1)) | ~isnan(offsetAB{2}(1,1))
handles = pushbutton17_Callback(hObject, eventdata, handles);
end
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
    filename488 = 'el_grII_corrected.ascii';
    elseif el532 == 1
    filename488 = 'el_grI_corrected.ascii';
    end
    FPName488 = [pathname filename488];
    filename642 = 'el_red_corrected.ascii';
    FPName642 = [pathname filename642];
    dlmwrite(FPName488, ABcorr{2});
    dlmwrite(FPName642, ABcorr{1});
end

handles.ABcorr = ABcorr;
guidata(hObject,handles);

catch errorObj
% If there is a problem, we display the error message
errordlg(getReport(errorObj,'extended','hyperlinks','off'),'Error');
end


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
try

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

catch errorObj
errordlg(getReport(errorObj,'extended','hyperlinks','off'),'Error');
end


% --- Executes on button press in pushbutton57.
function pushbutton57_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton57 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try

ABcorr = handles.ABcorr;
pixsize = str2double(get(handles.edit42, 'string'));

fov = FOV(ABcorr);
I = cell(1,2);
for i = 1:2
if ~isempty(ABcorr{i})
    I{i} = drpar(ABcorr{i}, pixsize, fov);
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
    filename488 = 'el_grII_corrected.ascii';
    elseif el532 == 1
    filename488 = 'el_grI_corrected.ascii';
    end
    FPName488 = [pathname filename488];
    filename642 = 'el_red_corrected.ascii';
    FPName642 = [pathname filename642];
    dlmwrite(FPName488, ABcorr{2});
    dlmwrite(FPName642, ABcorr{1});
end
handles.ABcorr = ABcorr;
handles.BW = BW;
guidata(hObject,handles);
end

catch errorObj
errordlg(getReport(errorObj,'extended','hyperlinks','off'),'Error');
end


% --- Executes on button press in pushbutton58.
function pushbutton58_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton58 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try

ABcorr = handles.ABcorr;
pixsize = str2double(get(handles.edit42, 'string'));


fov = FOV(ABcorr);
I = cell(2,1);
for i = 1:2
if (isempty(ABcorr{i}) == 0)
    I{i} = drpar(ABcorr{i}, pixsize, fov);
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
    filename488 = 'el_grII_corrected.ascii';
    elseif el532 == 1
    filename488 = 'el_grI_corrected.ascii';
    end
    FPName488 = [pathname filename488];
    filename642 = 'el_red_corrected.ascii';
    FPName642 = [pathname filename642];
    dlmwrite(FPName488, ABcorr{2});
    dlmwrite(FPName642, ABcorr{1});
end
handles.ABcorr = ABcorr;
guidata(hObject,handles);
end

catch errorObj
errordlg(getReport(errorObj,'extended','hyperlinks','off'),'Error');
end


% --- Executes on button press in pushbutton59.
function pushbutton59_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton59 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try
    
method = get(handles.popupmenu22, 'Value');
pix = str2double(get(handles.edit63, 'String'));
steps = str2double(get(handles.edit27, 'String'));
AB = handles.ABcorr;

FRC = cell(1,2);
resol = cell(1,2);
%parfor i = 1:2    
for i = 1:2  
if ~isempty (AB{i})
    [FRC{i}, resol{i}] = resolution(AB{i}, steps, method, pix);
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
           FPName = [pathname 'FRCred.png'];
        elseif el488 == 1
           FPName = [pathname 'FRCgrII.png'];
        elseif el532 == 1
           FPName = [pathname 'FRCgrI.png'];
        end
    h2 = figure('visible', 'off');
    if ~isempty (FRC{i})
    FRC1 = FRC{i};
    plot(FRC1(:,1),FRC1(:,3),'-',FRC1(:,1),FRC1(:,2),'.', 'LineWidth', 1, 'MarkerSize', 8);
    xlabel('Spatial frequency (nm^{-1})');
    ylabel('FRC');
    set(gca, 'XTick', 0:0.01:(0.5/pix), 'YTick', 0:0.2:1, 'FontSize', 14, 'LineWidth', 1);
    axis([0 (0.5/pix) 0 1]);
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

catch errorObj
errordlg(getReport(errorObj,'extended','hyperlinks','off'),'Error');
end

% --- Executes on button press in pushbutton60.
function pushbutton60_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton60 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try

FRC = handles.FRC;
pix = str2double(get(handles.edit63, 'String'));
if ~isempty (FRC{2})
    FRC1 = FRC{2};
    figure(2);
    plot(FRC1(:,1),FRC1(:,3),'-',FRC1(:,1),FRC1(:,2),'.', 'LineWidth', 1, 'MarkerSize', 8);
    xlabel('Spatial frequency, (nm^{-1})');
    ylabel('FRC');
    set(gca, 'XTick', 0:0.01:(0.5/pix), 'YTick', 0:0.2:1, 'FontSize', 14, 'LineWidth', 1);
    axis([0 (0.5/pix) 0 1]);
end

catch errorObj
errordlg(getReport(errorObj,'extended','hyperlinks','off'),'Error');
end


% --- Executes on button press in pushbutton62.
function pushbutton62_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton62 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try

FRC = handles.FRC;
pix = str2double(get(handles.edit63, 'String'));
if ~isempty (FRC{1})
    FRC1 = FRC{1};
    figure(1);
    plot(FRC1(:,1),FRC1(:,3),'-',FRC1(:,1),FRC1(:,2),'.', 'LineWidth', 1, 'MarkerSize', 8);
    xlabel('Spatial frequency, (nm^{-1})');
    ylabel('FRC');
    set(gca, 'XTick', 0:0.01:(0.5/pix), 'YTick', 0:0.2:1, 'FontSize', 14, 'LineWidth', 1);
    axis([0 (0.5/pix) 0 1]);
end

catch errorObj
errordlg(getReport(errorObj,'extended','hyperlinks','off'),'Error');
end


% --- Executes on button press in pushbutton63.
function pushbutton63_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton63 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try 
timeincolor(handles.ABcorr, handles.AB, handles.folder);

catch errorObj
errordlg(getReport(errorObj,'extended','hyperlinks','off'),'Error');
end



% --- Executes on button press in pushbutton77.
function pushbutton77_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton77 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try

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

catch errorObj
errordlg(getReport(errorObj,'extended','hyperlinks','off'),'Error');
end



% --- Executes on button press in pushbutton78.
function pushbutton78_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton78 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try

ABcorr = handles.ABcorr;
[~, ABcorr] = filtration(ABcorr);

save = get(handles.checkbox25, 'Value'); %auto save corrected eventlists
if save == 1
    el488 = get(handles.radiobutton11, 'Value');
    el532 = get(handles.radiobutton10, 'Value');
    pathname = handles.folder;
    if el488 == 1
    filename488 = 'el_grII_corrected.ascii';
    elseif el532 == 1
    filename488 = 'el_grI_corrected.ascii';
    end
    FPName488 = [pathname filename488];
    filename642 = 'el_red_corrected.ascii';
    FPName642 = [pathname filename642];
    dlmwrite(FPName488, ABcorr{2});
    dlmwrite(FPName642, ABcorr{1});
end

handles.ABcorr = ABcorr;
guidata(hObject,handles);

catch errorObj
errordlg(getReport(errorObj,'extended','hyperlinks','off'),'Error');
end


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
try 
Chroma_calibration;

catch errorObj
errordlg(getReport(errorObj,'extended','hyperlinks','off'),'Error');
end


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


% --- Executes on button press in checkbox58.
function checkbox58_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox58 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox58


% --- Executes on selection change in popupmenu22.
function popupmenu22_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu22 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu22


% --- Executes during object creation, after setting all properties.
function popupmenu22_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit63_Callback(hObject, eventdata, handles)
% hObject    handle to edit63 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit63 as text
%        str2double(get(hObject,'String')) returns contents of edit63 as a double


% --- Executes during object creation, after setting all properties.
function edit63_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit63 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu23.
function popupmenu23_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu23 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu23 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu23


% --- Executes during object creation, after setting all properties.
function popupmenu23_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu23 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function Plugins_Callback(hObject, eventdata, handles)
% hObject    handle to Plugins (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function clusters_Callback(hObject, eventdata, handles)
% hObject    handle to clusters (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isfield(handles, 'ABcorr')
    AB = handles.ABcorr;
    ClusterViSu(AB, handles.folder);
else
    AB = handles.AB;
    if ~isempty (AB{1}) || ~isempty (AB{2})
        ClusterViSu(AB, handles.folder);
    else
        ClusterViSu;
    end
end


% --------------------------------------------------------------------
function Untitled_1_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
