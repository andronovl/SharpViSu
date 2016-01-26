function varargout = Chroma_calibration(varargin)
% CHROMA_CALIBRATION MATLAB code for Chroma_calibration.fig
%      CHROMA_CALIBRATION, by itself, creates a new CHROMA_CALIBRATION or raises the existing
%      singleton*.
%
%      H = CHROMA_CALIBRATION returns the handle to a new CHROMA_CALIBRATION or the handle to
%      the existing singleton*.
%
%      CHROMA_CALIBRATION('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CHROMA_CALIBRATION.M with the given input arguments.
%
%      CHROMA_CALIBRATION('Property','Value',...) creates a new CHROMA_CALIBRATION or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Chroma_calibration_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Chroma_calibration_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Chroma_calibration

% Last Modified by GUIDE v2.5 01-Dec-2014 13:31:40

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Chroma_calibration_OpeningFcn, ...
                   'gui_OutputFcn',  @Chroma_calibration_OutputFcn, ...
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


% --- Executes just before Chroma_calibration is made visible.
function Chroma_calibration_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Chroma_calibration (see VARARGIN)

% Choose default command line output for Chroma_calibration
handles.output = hObject;
handles.ABC = cell(3,1);
handles.ABCcorr = cell(3,1);
handles.ABCcorr1 = cell(3,1);
handles.Anew = cell(2,1);
handles.Bnew = cell(2,1);

handles.loadfit = cell(2,1);
if isdeployed % Stand-alone mode.
    [~, result] = system('path');
    currentDir = char(regexpi(result, 'Path=(.*?);', 'tokens', 'once'));
else % MATLAB mode.
    currentDir = pwd;
end
if exist([currentDir '\532.mat'], 'file')
    handles.loadfit{1} = importdata([currentDir '\532.mat']);
    set(handles.text36, 'String', handles.loadfit{1}.Degree);
end
if exist([currentDir '\488.mat'], 'file')
    handles.loadfit{2} = importdata([currentDir '\488.mat']);
    set(handles.text37, 'String', handles.loadfit{2}.Degree);
end

handles.tform = cell(2,1);
handles.FOV = 18000;
if nargin == 1
    handles.folder = varargin{1};
else
    handles.folder = currentDir;
end
% Update handles structure
guidata(hObject, handles);
if ~isempty(handles.loadfit{1}) || ~isempty(handles.loadfit{2})
    UpdateImage(hObject, handles);
end

% UIWAIT makes Chroma_calibration wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Chroma_calibration_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
GetFile(hObject, handles, 1);
handles.FOV = FOV(handles.ABC);
set(handles.text14, 'String', handles.FOV/1000);


function GetFile( hObject, handles, channel)
format = get(handles.popupmenu1, 'Value');
    extension = '*.*';
if isfield(handles, 'folder')
    pathname = handles.folder;
    [filename, pathname] = uigetfile({extension}, 'Please select the eventlist of an aquisition with beads', pathname);
else
    [filename, pathname] = uigetfile({extension}, 'Please select the eventlist of an aquisition with beads');
end
if pathname ~= 0
    if channel == 1 % 642
        set (handles.edit2, 'String', [pathname, filename]);
    elseif channel == 2 % 532
        set (handles.edit4, 'String', [pathname, filename]);
    elseif channel == 3 % 488
        set (handles.edit5, 'String', [pathname, filename]);
    end
    
A = importdata([pathname, filename]);
if isstruct(A)
    A = A.data;
end
A = adjustformat(A, format);    % adjust the format to the standard

handles.folder = pathname;
handles.ABC{channel} = A;
else
    if channel == 1 % 642
        set (handles.edit2, 'String', []);
        handles.Anew{1} = [];
        handles.Anew{2} = [];
    elseif channel == 2 % 532
        set (handles.edit4, 'String', []);
        handles.Bnew{1} = [];
    elseif channel == 3 % 488
        set (handles.edit5, 'String', []);
        handles.Bnew{2} = [];
    end
    handles.ABC{channel} = [];
end
handles.ABCcorr{channel} = handles.ABC{channel};
handles.ABCcorr1{channel} = handles.ABC{channel};
handles.FOV = FOV(handles.ABC);
guidata(hObject,handles);
updatehist(hObject, handles);



% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
GetFile(hObject, handles, 2);
handles.FOV = FOV(handles.ABC);
set(handles.text14, 'String', handles.FOV/1000);



% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
GetFile(hObject, handles, 3);
handles.FOV = FOV(handles.ABC);
set(handles.text14, 'String', handles.FOV/1000);



function edit10_Callback(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Fitting(hObject, handles);



% --- Executes on button press in pushbutton7.
function pushbutton7_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
tform = handles.tform;
%function currentDir = getcurrentdir http://www.mathworks.com/matlabcentral/answers/92949-how-can-i-find-the-directory-containing-my-compiled-application
if isdeployed % Stand-alone mode.
    [~, result] = system('path');
    currentDir = char(regexpi(result, 'Path=(.*?);', 'tokens', 'once'));
else % MATLAB mode.
    currentDir = pwd;
end
%dir = uigetdir(currentDir, 'Installation folder of the program');
%dir = [dir '\'];
%if FileName ~= 0
FileName{1} = '\532.mat';
FileName{2} = '\488.mat';
for i = 1 : 2
    if ~isempty(tform{i})
        a = tform{i};
FPName=[currentDir FileName{i}];
save(FPName, 'a');
    end
end
if isempty(tform{1}) && isempty(tform{2})
    h = errordlg('Please load calibration data first','No calibration data');
end



function edit14_Callback(hObject, eventdata, handles)
% hObject    handle to edit14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
filter_photons(hObject,handles);

function filter_photons(hObject, handles)
try
    
ABC = handles.ABCcorr;

beg(1) = 1000 * str2double(get(handles.edit14, 'String')); % lowest photons A
fin(1) = Inf; % highest photons B

beg(2) = 1000 * str2double(get(handles.edit15, 'String')); % lowest photons B
fin(2) = Inf; % highest photons B

beg(3) = 1000 * str2double(get(handles.edit16, 'String')); % lowest photons B
fin(3) = Inf; % highest photons B


for i = 1:3
if ~isempty(ABC{i})
  ABC{i} = keepphotons(ABC{i}, beg(i), fin(i));
end
end
handles.ABCcorr1 = ABC;
guidata(hObject,handles);
updatehist(hObject, handles);
Fitting(hObject, handles);

catch errorObj
    errordlg(getReport(errorObj,'extended','hyperlinks','off'),'Error');
end



function edit15_Callback(hObject, eventdata, handles)
% hObject    handle to edit15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
filter_photons(hObject,handles);


function edit16_Callback(hObject, eventdata, handles)
% hObject    handle to edit16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
filter_photons(hObject,handles);



% --- Executes on selection change in popupmenu2.
function popupmenu2_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
UpdatePlot(hObject,handles);



% --- Executes on button press in checkbox3.
function checkbox3_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
UpdatePlot(hObject,handles);



function edit17_Callback(hObject, eventdata, handles)
% hObject    handle to edit17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
updatehist(hObject, handles);


function edit18_Callback(hObject, eventdata, handles)
% hObject    handle to edit18 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
UpdateImage(hObject, handles);


% --- Executes on selection change in popupmenu4.
function popupmenu4_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
UpdateImage(hObject, handles);


% --- Executes during object creation, after setting all properties.
function popupmenu4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu5.
function popupmenu5_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
UpdateImage(hObject, handles);


% --- Executes during object creation, after setting all properties.
function popupmenu5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit19_Callback(hObject, eventdata, handles)
% hObject    handle to edit19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
UpdateImage(hObject, handles);



% --- Executes on button press in pushbutton9.
function pushbutton9_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try
rad(1) = str2double(get(handles.edit7, 'string')); % radius A
gap(1) = str2double(get(handles.edit20, 'string')); % nb of empty frames
rad(2) = str2double(get(handles.edit8, 'string')); % radius B
gap(2) = str2double(get(handles.edit21, 'string')); % nb of empty frames
rad(3) = str2double(get(handles.edit9, 'string')); % radius C
gap(3) = str2double(get(handles.edit22, 'string')); % nb of empty frames

ABC = handles.ABC;
ABCcorr = cell(3,1);
for i = 1:3
if ~isempty(ABC{i}) && rad(i) ~= 0
  ABCcorr{i} = filtercons(ABC{i}, rad(i), gap(i), 0);
else
  ABCcorr{i} = ABC{i};
end
end
if isempty(ABC{1}) && isempty(ABC{2}) && isempty(ABC{3})
    h = errordlg('Please load calibration data first','No calibration data');
end

handles.ABCcorr = ABCcorr;
handles.ABCcorr1 = ABCcorr;
guidata(hObject,handles);

filter_photons(hObject,handles);

catch errorObj
    errordlg(getReport(errorObj,'extended','hyperlinks','off'),'Error');
end


function updatehist(hObject, handles)
ABC = handles.ABCcorr1;
show = get(handles.popupmenu6, 'Value');
if ~isempty(ABC{1})
    la = size(ABC{1},1);
    set(handles.text19, 'String', la);
if show == 1
    bins = str2double(get(handles.edit17, 'String')); % nb of bins for red histogram
    hist(handles.axes1, ABC{1}(:,7), bins);
end
end

if ~isempty(ABC{2})
    la = size(ABC{2},1);
    set(handles.text20, 'String', la);
    if show == 2
    bins = str2double(get(handles.edit17, 'String')); % nb of bins for red histogram
    hist(handles.axes1, ABC{2}(:,7), bins);
    end
end
if ~isempty(ABC{3})
    la = size(ABC{3},1);
    set(handles.text21, 'String', la);
    
if show == 3
    bins = str2double(get(handles.edit17, 'String')); % nb of bins for red histogram
    hist(handles.axes1, ABC{3}(:,7), bins);
end
end
UpdateImage(hObject, handles)



function A = keepphotons (A, beg, fin)
% keeps only events wigh photon count from beg to fin
A(((A(:,7) < beg) | (A(:,7) > fin)),:) = [];


% --- Executes on selection change in popupmenu6.
function popupmenu6_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
updatehist(hObject, handles)



function edit23_Callback(hObject, eventdata, handles)
% hObject    handle to edit23 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
UpdateImage(hObject, handles);


% --- Executes when selected object is changed in uipanel4.
function uipanel4_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel4 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
UpdateImage(hObject, handles);



function UpdateImage(hObject, handles)
try
ABC = handles.ABC;
ABCcorr = handles.ABCcorr1;
Anew = handles.Anew;
Bnew = handles.Bnew;
show = get(handles.popupmenu7, 'Value');
fov = handles.FOV/1000;
if get(handles.radiobutton1, 'Value') % if show scatter
if (isempty(Anew{1}) && isempty(Anew{2})) || show == 1 || show == 2
    cla(handles.axes2);
    if show == 1 % show initial events
    hold on
    scatter(ABC{1}(:,4)/1000, ABC{1}(:,5)/1000, '.r', 'Parent', handles.axes2);
    if ~isempty(ABC{2})
    scatter(ABC{2}(:,4)/1000, ABC{2}(:,5)/1000, '.g', 'Parent', handles.axes2);
    end
    if ~isempty(ABC{3})
    scatter(ABC{3}(:,4)/1000, ABC{3}(:,5)/1000, '.b', 'Parent', handles.axes2);
    end
    hold off
    else % show filtered events
    hold on
    scatter(ABCcorr{1}(:,4)/1000, ABCcorr{1}(:,5)/1000, '.r', 'Parent', handles.axes2);
    if ~isempty(ABCcorr{2})
    scatter(ABCcorr{2}(:,4)/1000, ABCcorr{2}(:,5)/1000, '.g', 'Parent', handles.axes2);
    end
    if ~isempty(ABCcorr{3})
    scatter(ABCcorr{3}(:,4)/1000, ABCcorr{3}(:,5)/1000, '.b', 'Parent', handles.axes2);
    end
    hold off
    end
end
if ~isempty(Anew{1}) || ~isempty(Anew{2})
    if show == 3 % show events used for the fit
    cla(handles.axes2);
    hold on
    if ~isempty (Anew{1})
    scatter(Anew{1}(:,1)/1000, Anew{1}(:,2)/1000, 20, 'sr', 'Parent', handles.axes2, 'MarkerFaceColor', 'r');
    scatter(Bnew{1}(:,1)/1000, Bnew{1}(:,2)/1000, 20, 'sg', 'Parent', handles.axes2, 'MarkerFaceColor', 'g');
    end
    if ~isempty (Anew{2})
    scatter(Anew{2}(:,1)/1000, Anew{2}(:,2)/1000, 30, 'dr', 'Parent', handles.axes2, 'MarkerFaceColor', 'r');
    scatter(Bnew{2}(:,1)/1000, Bnew{2}(:,2)/1000, 30, 'db', 'Parent', handles.axes2, 'MarkerFaceColor', 'b');
    end
    hold off
    end
    if show == 4 % show the fit
        if ~get(handles.checkbox4, 'Value')
        cla(handles.axes2);
        tform = handles.tform;
    hold on
    if ~isempty (Anew{1})
    Bcorr{1} = transformPointsInverse(tform{1}, Bnew{1}/1000);
    scatter(Anew{1}(:,1)/1000, Anew{1}(:,2)/1000, 20, 'sr', 'Parent', handles.axes2, 'MarkerFaceColor', 'r');
    scatter(Bcorr{1}(:,1), Bcorr{1}(:,2), 20, 'sg', 'Parent', handles.axes2, 'MarkerFaceColor', 'g');
    end
    if ~isempty (Anew{2})
    Bcorr{2} = transformPointsInverse(tform{2}, Bnew{2}/1000);
    scatter(Anew{2}(:,1)/1000, Anew{2}(:,2)/1000, 30, 'dr', 'Parent', handles.axes2, 'MarkerFaceColor', 'r');
    scatter(Bcorr{2}(:,1), Bcorr{2}(:,2), 30, 'db', 'Parent', handles.axes2, 'MarkerFaceColor', 'b');
    end
    hold off
        end
        if get(handles.checkbox4, 'Value') % if use a calibration
            cla(handles.axes2);
        tform = handles.loadfit;
        hold on
    if ~isempty (Anew{1}) && ~isempty (tform{1})
    Bcorr{1} = transformPointsInverse(tform{1}, Bnew{1}/1000);
    scatter(Anew{1}(:,1)/1000, Anew{1}(:,2)/1000, 20, 'sr', 'Parent', handles.axes2, 'MarkerFaceColor', 'r');
    scatter(Bcorr{1}(:,1), Bcorr{1}(:,2), 20, 'sg', 'Parent', handles.axes2, 'MarkerFaceColor', 'g');
    end
    if ~isempty (Anew{2}) && ~isempty (tform{2})
    Bcorr{2} = transformPointsInverse(tform{2}, Bnew{2}/1000);
    scatter(Anew{2}(:,1)/1000, Anew{2}(:,2)/1000, 30, 'dr', 'Parent', handles.axes2, 'MarkerFaceColor', 'r');
    scatter(Bcorr{2}(:,1), Bcorr{2}(:,2), 30, 'db', 'Parent', handles.axes2, 'MarkerFaceColor', 'b');
    end
    hold off
        end
    end
end
% % save the scatter plot
% 
% pathname = handles.folder;
% FPName = [pathname 'scatter_for fitting.png'];
% h2 = figure('visible', 'off');
% 
% 
% 
% hold on
% % if ~isempty (Anew{1})
% % scatter(Anew{1}(:,1)/1000, Anew{1}(:,2)/1000, 50, 'xr', 'MarkerFaceColor', 'r');
% % scatter(Bnew{1}(:,1)/1000, Bnew{1}(:,2)/1000, 50, 'xg', 'MarkerFaceColor', 'g');
% % end
% if ~isempty (Anew{2})
% scatter(Anew{2}(:,1)/1000, Anew{2}(:,2)/1000, 500, '+r', 'MarkerFaceColor', 'r', 'LineWidth', 1.5);
% scatter(Bnew{2}(:,1)/1000, Bnew{2}(:,2)/1000, 500, 'xg', 'MarkerFaceColor', 'g', 'LineWidth', 1.5);
% end
% hold off
% 
% 
% xlabel('X, µm', 'FontSize', 18);
% ylabel('Y, µm', 'FontSize', 18);
% set(gca, 'xaxislocation', 'top', 'ydir', 'reverse', 'FontSize', 18, 'LineWidth', 0.5, 'DataAspectRatio', [1 1 1], 'XLim', [0 fov], 'YLim', [0 fov], 'XTick', [0 6 12 18], 'YTick', [0 6 12 18]);
% %axis(gca, [0 fov 0 fov]);
% box on;
% myStyle = hgexport('factorystyle');
% myStyle.Format = 'png';
% myStyle.Width = 6;
% myStyle.Height = 6;
% myStyle.Resolution = 300;
% myStyle.Units = 'inch';
% myStyle.FixedFontSize = 8;
% hgexport(h2, FPName, myStyle, 'Format', 'png');
% close(h2);

elseif get(handles.radiobutton2, 'Value') % if show vectors

cla(handles.axes2);
if get(handles.checkbox4, 'Value') % if use a calibration
tform = handles.loadfit;
else
tform = handles.tform;
end
step = str2double(get(handles.edit19, 'String'));
length = str2double(get(handles.edit25, 'String'));
origin = get(handles.popupmenu4, 'Value');
tip = get(handles.popupmenu5, 'Value');
[X,Y] = meshgrid(0:step:fov);
X = reshape(X, numel(X), 1);
Y = reshape(Y, numel(Y), 1);
XY = [X, Y];
if ~isempty (tform{1}) && (origin == 2 || tip == 2) % 532 corr
XYcorr532 = transformPointsInverse(tform{1}, XY); 
end
if ~isempty (tform{2}) && (origin == 3 || tip == 3) % 488 corr
XYcorr488 = transformPointsInverse(tform{2}, XY); 
end

if origin == 1 && tip == 2 % 642-532
    XYshift = (XY - XYcorr532) * length;
    quiver('Parent', handles.axes2, XYcorr532(:,1), XYcorr532(:,2), XYshift(:,1), XYshift(:,2), 0, 'LineWidth', 1.5, 'MaxHeadSize', 0.5, 'Color', [0.5 0.5 0]);
elseif origin == 1 && tip == 3 % 642-488
    XYshift = (XY - XYcorr488) * length;
    quiver('Parent', handles.axes2, XYcorr488(:,1), XYcorr488(:,2), XYshift(:,1), XYshift(:,2), 0, 'LineWidth', 1.5, 'MaxHeadSize', 0.5, 'Color', [0.5 0 0.5]);
elseif origin == 2 && tip == 1 % 532-642
    XYshift = (XYcorr532 - XY) * length;
    quiver('Parent', handles.axes2, XY(:,1), XY(:,2), XYshift(:,1), XYshift(:,2), 0, 'LineWidth', 1.5, 'MaxHeadSize', 0.5, 'Color', [0.5 0.5 0]);
elseif origin == 2 && tip == 3 % 532-488
    XYshift = (XYcorr532 - XYcorr488) * length;
    quiver('Parent', handles.axes2, XYcorr532(:,1), XYcorr532(:,2), XYshift(:,1), XYshift(:,2), 0, 'LineWidth', 1.5, 'MaxHeadSize', 0.5, 'Color', [0 0.5 0.5]);
elseif origin == 3 && tip == 1 % 488-642
    XYshift = (XYcorr488 - XY) * length;
    quiver('Parent', handles.axes2, XY(:,1), XY(:,2), XYshift(:,1), XYshift(:,2), 0, 'LineWidth', 1.5, 'MaxHeadSize', 0.5, 'Color', [0.5 0 0.5]);
elseif origin == 3 && tip == 2 % 488-532
    XYshift = (XYcorr488 - XYcorr532) * length;
    quiver('Parent', handles.axes2, XYcorr488(:,1), XYcorr488(:,2), XYshift(:,1), XYshift(:,2), 0, 'LineWidth', 1.5, 'MaxHeadSize', 0.5, 'Color', [0 0.5 0.5]);

end
    
end
set(handles.axes2,'ydir','reverse', 'color', [0.9 0.9 0.9], 'FontSize', 10);
axis(handles.axes2, [0 fov 0 fov]);

% %save the quiver plot
% pathname = handles.folder;
% FPName = [pathname 'vectors_642_488.png'];
% h2 = figure('visible', 'off');
% quiver(XYcorr488(:,1), XYcorr488(:,2), XYshift(:,1), XYshift(:,2), 0, 'LineWidth', 1.5, 'MaxHeadSize', 0.5, 'Color', [0.5 0 0.5]);
% xlabel('X, µm', 'FontSize', 18);
% ylabel('Y, µm', 'FontSize', 18);
% set(gca, 'xaxislocation', 'top', 'ydir', 'reverse', 'FontSize', 18, 'LineWidth', 0.5, 'DataAspectRatio', [1 1 1], 'XLim', [0 fov], 'YLim', [0 fov], 'XTick', [0 6 12 18], 'YTick', [0 6 12 18]);
% %axis(gca, [0 fov 0 fov]);
% myStyle = hgexport('factorystyle');
% myStyle.Format = 'png';
% myStyle.Width = 6;
% myStyle.Height = 6;
% myStyle.Resolution = 300;
% myStyle.Units = 'inch';
% myStyle.FixedFontSize = 8;
% hgexport(h2, FPName, myStyle, 'Format', 'png');
% close(h2);
guidata(hObject,handles);

catch errorObj
    errordlg(getReport(errorObj,'extended','hyperlinks','off'),'Error');
end


function edit24_Callback(hObject, eventdata, handles)
% hObject    handle to edit24 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Fitting(hObject, handles);



function Fitting(hObject, handles)
try
    
ABC = handles.ABCcorr;

R = str2double(get(handles.edit24, 'String')); % search radius
order = str2double(get(handles.edit10, 'String')); % search radius
tform = cell(2,1);
Anew = cell(2,1);
Bnew = cell(2,1);
for i = 1:2
    if ~isempty(ABC{i+1})
dist = pdist2(ABC{1}(:,4:5), ABC{i+1}(:,4:5));

Bin = dist < R; % find closets pairs of points
% delete points if there are more than 1 point from the other channel within R
S1 = find(sum(Bin,1) > 1); % index of the column to put all zeros
S2 = find(sum(Bin,2) > 1); % index of the row to put all zeros
if ~isempty(S2)
Bin(S2,:) = zeros(size(S2,1),size(Bin,2));
end
if ~isempty(S1)
Bin(:,S1) = zeros(size(Bin,1),size(S1,2));
end
[a,b] = ind2sub(size(Bin),find(Bin)); %indices of A, indices of B(C)
Anew{i} = ABC{1}(a,4:5);
Bnew{i} = ABC{i+1}(b,4:5);
tform{i} = fitgeotrans(Anew{i}/1000,Bnew{i}/1000,'polynomial',order);
    end
end

handles.tform = tform;
handles.Anew = Anew;
handles.Bnew = Bnew;
guidata(hObject,handles);
UpdateImage(hObject, handles);

catch errorObj
    errordlg(getReport(errorObj,'extended','hyperlinks','off'),'Error');
end



% --- Executes on selection change in popupmenu7.
function popupmenu7_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
UpdateImage(hObject, handles);


% --- Executes on button press in pushbutton10.
function pushbutton10_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
extension = '*.mat*';
if isfield(handles, 'folder')
    pathname = handles.folder;
    [filename, pathname] = uigetfile({extension}, 'Please select a calibration file for 532 nm channel', pathname);
else
    [filename, pathname] = uigetfile({extension}, 'Please select a calibration file for 532 nm channel');
end
if pathname ~= 0
    
handles.loadfit{1} = importdata([pathname, filename]);
set(handles.text36, 'String', handles.loadfit{1}.Degree);

end

if pathname ~= 0
    [filename, pathname] = uigetfile({extension}, 'Please select a calibration file for 488 nm channel', pathname);
else
    [filename, pathname] = uigetfile({extension}, 'Please select a calibration file for 488 nm channel');
end
if pathname ~= 0
    
handles.loadfit{2} = importdata([pathname, filename]);
set(handles.text37, 'String', handles.loadfit{2}.Degree);

end
guidata(hObject, handles);
UpdateImage(hObject, handles);


% --- Executes on button press in checkbox4.
function checkbox4_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
UpdateImage(hObject, handles);



function edit25_Callback(hObject, eventdata, handles)
% hObject    handle to edit25 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
UpdateImage(hObject, handles);
