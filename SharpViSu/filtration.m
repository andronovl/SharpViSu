function varargout = filtration(varargin)
% FILTRATION MATLAB code for filtration.fig
%      FILTRATION, by itself, creates a new FILTRATION or raises the existing
%      singleton*.
%
%      H = FILTRATION returns the handle to a new FILTRATION or the handle to
%      the existing singleton*.
%
%      FILTRATION('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FILTRATION.M with the given input arguments.
%
%      FILTRATION('Property','Value',...) creates a new FILTRATION or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before filtration_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to filtration_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help filtration

% Last Modified by GUIDE v2.5 23-Dec-2014 12:04:59

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @filtration_OpeningFcn, ...
                   'gui_OutputFcn',  @filtration_OutputFcn, ...
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

% --- Executes just before filtration is made visible.
function filtration_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to filtration (see VARARGIN)

% Choose default command line output for filtration
handles.output = hObject;
AB = varargin{1};
ABcorr = AB;
handles.AB = AB;
handles.ABcorr = ABcorr;

% evaluate initial eventlist
updatefigure(handles, ABcorr);

% Update handles structure
guidata(hObject, handles);

% This sets up the initial plot - only do when we are invisible
% so window can get raised using filtration.

% UIWAIT makes filtration wait for user response (see UIRESUME)
uiwait(hObject);


% --- Outputs from this function are returned to the command line.
function varargout = filtration_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
if isfield(handles, 'ABcorr')
varargout{2} = handles.ABcorr;
else
varargout{2} = handles.AB;
end
delete(hObject);



% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try
if ~get(handles.checkbox3, 'Value')
    AB = handles.AB; % correct exclusively initial eventlists
else
    AB = handles.ABcorr; % correct already filtered eventlists
end

rad(1) = str2double(get(handles.edit1, 'string')); % radius A
gap(1) = str2double(get(handles.edit2, 'string')); % nb of empty frames
anis(1) = get(handles.checkbox1, 'Value'); % anisotropically
rad(2) = str2double(get(handles.edit3, 'string')); % radius B
gap(2) = str2double(get(handles.edit4, 'string')); % nb of empty frames
anis(2) = get(handles.checkbox2, 'Value'); % anisotropically

ABcorr = handles.ABcorr;
for i = 1:2
if ~isempty(AB{i}) && rad(i) ~= 0
  ABcorr{i} = filtercons(AB{i}, rad(i), gap(i), anis(i));
else
  ABcorr{i} = AB{i};
end
end
updatefigure(handles, ABcorr);

handles.ABcorr = ABcorr;
guidata(hObject,handles);

catch errorObj
    errordlg(getReport(errorObj,'extended','hyperlinks','off'),'Error');
end





function edit9_Callback(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of edit9 as text
%        str2double(get(hObject,'String')) returns contents of edit9 as a double
AB = handles.ABcorr{1};

beg = str2double(get(handles.edit9, 'String')); % lowest photons
fin = str2double(get(handles.edit10, 'String')); % highest photons

if ~isempty(AB)
  AB = keepphotons(AB, beg, fin);
end
ABcell = cell(1,2);
ABcell{1} = AB;
updatefigure(handles, ABcell);





function edit10_Callback(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of edit10 as text
%        str2double(get(hObject,'String')) returns contents of edit10 as a double
AB = handles.ABcorr{1};

beg = str2double(get(handles.edit9, 'String')); % lowest photons
fin = str2double(get(handles.edit10, 'String')); % highest photons

if ~isempty(AB)
  AB = keepphotons(AB, beg, fin);
end
ABcell = cell(1,2);
ABcell{1} = AB;
updatefigure(handles, ABcell);




function edit11_Callback(hObject, eventdata, handles)
% hObject    handle to edit11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of edit11 as text
%        str2double(get(hObject,'String')) returns contents of edit11 as a double
AB = handles.ABcorr{2};

beg = str2double(get(handles.edit11, 'String')); % lowest photons
fin = str2double(get(handles.edit12, 'String')); % highest photons

if ~isempty(AB)
  AB = keepphotons(AB, beg, fin);
end
ABcell = cell(1,2);
ABcell{2} = AB;
updatefigure(handles, ABcell);




function edit12_Callback(hObject, eventdata, handles)
% hObject    handle to edit12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of edit12 as text
%        str2double(get(hObject,'String')) returns contents of edit12 as a double
AB = handles.ABcorr{2};

beg = str2double(get(handles.edit11, 'String')); % lowest photons
fin = str2double(get(handles.edit12, 'String')); % highest photons

if ~isempty(AB)
  AB = keepphotons(AB, beg, fin);
end
ABcell = cell(1,2);
ABcell{2} = AB;
updatefigure(handles, ABcell);




% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
AB = handles.ABcorr;

beg(1) = str2double(get(handles.edit9, 'String')); % lowest photons A
fin(1) = str2double(get(handles.edit10, 'String')); % highest photons B

beg(2) = str2double(get(handles.edit11, 'String')); % lowest photons B
fin(2) = str2double(get(handles.edit12, 'String')); % highest photons B

for i = 1:2
if ~isempty(AB{i})
  AB{i} = keepphotons(AB{i}, beg(i), fin(i));
end
end
updatefigure(handles, AB);

handles.ABcorr = AB;
guidata(hObject,handles);



function edit5_Callback(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of edit5 as text
%        str2double(get(hObject,'String')) returns contents of edit5 as a double
AB = handles.ABcorr{1};

beg = str2double(get(handles.edit5, 'string')); % first frame
fin = str2double(get(handles.edit6, 'string')); % last frame

if ~isempty(AB)
  AB = keepframes(AB, beg, fin);
end
ABcell = cell(1,2);
ABcell{1} = AB;
updatefigure(handles, ABcell);




function edit6_Callback(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of edit6 as text
%        str2double(get(hObject,'String')) returns contents of edit6 as a double
AB = handles.ABcorr{1};

beg = str2double(get(handles.edit5, 'string')); % first frame
fin = str2double(get(handles.edit6, 'string')); % last frame

if ~isempty(AB)
  AB = keepframes(AB, beg, fin);
end
ABcell = cell(1,2);
ABcell{1} = AB;
updatefigure(handles, ABcell);





function edit7_Callback(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of edit7 as text
%        str2double(get(hObject,'String')) returns contents of edit7 as a double
AB = handles.ABcorr{2};

beg = str2double(get(handles.edit7, 'string')); % first frame
fin = str2double(get(handles.edit8, 'string')); % last frame

if ~isempty(AB)
  AB = keepframes(AB, beg, fin);
end
ABcell = cell(1,2);
ABcell{2} = AB;
updatefigure(handles, ABcell);





function edit8_Callback(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of edit8 as text
%        str2double(get(hObject,'String')) returns contents of edit8 as a double
AB = handles.ABcorr{2};

beg = str2double(get(handles.edit7, 'string')); % first frame
fin = str2double(get(handles.edit8, 'string')); % last frame

if ~isempty(AB)
  AB = keepframes(AB, beg, fin);
end
ABcell = cell(1,2);
ABcell{2} = AB;
updatefigure(handles, ABcell);




% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
AB = handles.ABcorr;

beg(1) = str2double(get(handles.edit5, 'string')); % first frame A
fin(1) = str2double(get(handles.edit6, 'string')); % last frame B

beg(2) = str2double(get(handles.edit7, 'string')); % first frame B
fin(2) = str2double(get(handles.edit8, 'string')); % last frame B

for i = 1:2
if ~isempty(AB{i})
  AB{i} = keepframes(AB{i}, beg(i), fin(i));
end
end
updatefigure(handles, AB);

handles.ABcorr = AB;
guidata(hObject,handles);




% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
if isequal(get(hObject, 'waitstatus'), 'waiting')
    % The GUI is still in UIWAIT, use UIRESUME
    uiresume(hObject);
else
    % The GUI is no longer waiting, just close it
    delete(hObject);
end



function plotevperfr(A, h)
%plots function events per frame
% A eventlist, h handle to axes
evperfr = eventperframe(A);
plot(h, evperfr);
set(h, 'XMinorTick', 'on', 'FontSize', 9);
title(h, 'Events per frame', 'FontSize', 9);
xlabel(h, 'Frame number', 'FontSize', 9);
ylabel(h, 'Number of events', 'FontSize', 9);


function plotphperev(A, bins, h, hmean, hmedian, hmode)
%plots histogram of the distibution number of photones in a frame
% A eventlist, h handle to axes, b number of bins of hisogram
hist(h, A(:,7), bins);
set(h, 'XMinorTick', 'on', 'FontSize', 9);
xlabel(h, 'Number of photons', 'FontSize', 9);
ylabel(h, 'Number of events', 'FontSize', 9);

if exist ('hmean', 'var')
set(hmean, 'String', mean(A(:,7)));
set(hmedian, 'String', median(A(:,7)));
set(hmode, 'String', mode(A(:,7)));
end
% save the histogram

% pathname = 'D:\My Documents\MATLAB\20130301_tub_al647\';
% FPName = [pathname 'histPhotons.png'];
% h2 = figure('visible', 'off');
% hist(A(:,7), bins);
% %n=get(gca,'xtick');
% %set(gca,'xticklabel',sprintf('%2.0e|',n'));
% xlabel('Number of photons', 'FontSize', 18);
% ylabel('Number of localizations', 'FontSize', 18);
% %set(gca, 'FontSize', 18, 'LineWidth', 0.5, 'DataAspectRatio', [1 1 1], 'XLim', [0 fov], 'YLim', [0 fov], 'XTick', [0 6 12 18], 'YTick', [0 6 12 18]);
% set(gca, 'FontSize', 18, 'LineWidth', 0.5);
% %axis(gca, [0 fov 0 fov]);
% myStyle = hgexport('factorystyle');
% myStyle.Format = 'png';
% myStyle.Width = 10;
% myStyle.Height = 6;
% myStyle.Resolution = 300;
% myStyle.Units = 'inch';
% myStyle.FixedFontSize = 8;
% hgexport(h2, FPName, myStyle, 'Format', 'png');
% close(h2);


function edit13_Callback(hObject, eventdata, handles)
% hObject    handle to edit13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of edit13 as text
%        str2double(get(hObject,'String')) returns contents of edit13 as a double
AB = handles.ABcorr;
if ~isempty(AB{1})
    bins = str2double(get(hObject, 'String')); % nb of bins for red histogram
    plotphperev(AB{1}, bins, handles.axes1)
end


function edit14_Callback(hObject, eventdata, handles)
% hObject    handle to edit14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of edit14 as text
%        str2double(get(hObject,'String')) returns contents of edit14 as a double
AB = handles.ABcorr;
if ~isempty(AB{2})
    bins = str2double(get(hObject, 'String')); % nb of bins for red histogram
    plotphperev(AB{2}, bins, handles.axes2)
end



function updatefigure(handles, ABcorr)
% ABcorr = cell 
try
if ~isempty(ABcorr{1})
    plotevperfr(ABcorr{1}, handles.axes5);
    bins = str2double(get(handles.edit13, 'String')); % nb of bins for red histogram
    plotphperev(ABcorr{1}, bins, handles.axes1, handles.text25, handles.text26, handles.text27)
    set(handles.text51, 'String', size(ABcorr{1}, 1));
end
if ~isempty(ABcorr{2})
    plotevperfr(ABcorr{2}, handles.axes6);
    bins = str2double(get(handles.edit14, 'String')); % nb of bins for green histogram
    plotphperev(ABcorr{2}, bins, handles.axes2, handles.text46, handles.text47, handles.text48)
    set(handles.text52, 'String', size(ABcorr{2}, 1));
end

catch errorObj
    errordlg(getReport(errorObj,'extended','hyperlinks','off'),'Error');
end
        
function A = keepframes (A, beg, fin)
% keeps only events from frames from beg to fin
try
A(((A(:,2) < beg) | (A(:,2) > fin)),:) = [];

catch errorObj
    errordlg(getReport(errorObj,'extended','hyperlinks','off'),'Error');
end

function A = keepphotons (A, beg, fin)
% keeps only events wigh photon count from beg to fin
try
A(((A(:,7) < beg) | (A(:,7) > fin)),:) = [];

catch errorObj
    errordlg(getReport(errorObj,'extended','hyperlinks','off'),'Error');
end
