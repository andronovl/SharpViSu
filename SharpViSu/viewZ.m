function varargout = viewZ(varargin)
% VIEWZ MATLAB code for viewZ.fig
%      VIEWZ, by itself, creates a new VIEWZ or raises the existing
%      singleton*.
%
%      H = VIEWZ returns the handle to a new VIEWZ or the handle to
%      the existing singleton*.
%
%      VIEWZ('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in VIEWZ.M with the given input arguments.
%
%      VIEWZ('Property','Value',...) creates a new VIEWZ or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before viewZ_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to viewZ_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help viewZ

% Last Modified by GUIDE v2.5 14-Jan-2015 14:24:40

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @viewZ_OpeningFcn, ...
                   'gui_OutputFcn',  @viewZ_OutputFcn, ...
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

% --- Executes just before viewZ is made visible.
function viewZ_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to viewZ (see VARARGIN)

% Choose default command line output for viewZ
handles.output = hObject;

AB = varargin{1};
handles.AB = AB;
guidata(hObject, handles);
GetSlices(hObject, handles);

addlistener(handles.slider1,'ContinuousValueChange',@slider_Callback);

% UIWAIT makes viewZ wait for user response (see UIRESUME)
uiwait(handles.figure);


% --- Outputs from this function are returned to the command line.
function varargout = viewZ_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
varargout{2} = handles.RGB;
%varargout{3} = handles.multiplier;
delete(handles.figure);




% --- Executes when user attempts to close figure.
function figure_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
if isequal(get(handles.figure, 'waitstatus'), 'waiting')
    % The GUI is still in UIWAIT, use UIRESUME
    uiresume(handles.figure);
else
    % The GUI is no longer waiting, just close it
     delete(handles.figure);
end



function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double
updateimage(hObject, handles);



function slider_Callback(hObject, eventdata, handles)
handles = guidata(hObject);
slice = get(handles.slider1, 'Value');
set(handles.edit3, 'String', round(slice));
guidata(hObject, handles);
ShowSlice(hObject, handles);




function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double
slider_value = str2double(get(handles.edit3,'String'));
set(handles.slider1, 'Value', slider_value);
ShowSlice(hObject, handles);





function updateimage(hObject, handles)
% updates the image
try
I0 = handles.I0;
mult(1) = str2double(get(handles.edit6, 'String'));
mult(2) = str2double(get(handles.edit7, 'String'));
value = get(handles.popupmenu6, 'Value');
stringlist = get(handles.popupmenu6, 'String');
mapstr{1} = stringlist{value};
value = get(handles.popupmenu7, 'Value');
stringlist = get(handles.popupmenu7, 'String');
mapstr{2} = stringlist{value};

I = cell(size(I0));
RGB0 = cell(size(I0));
map = cell(2,1);
for i = 1 : 2
    I(i,:) = cellfun(@(x) x * mult(i), I0(i,:), 'UniformOutput', 0);
    if strcmp(mapstr{i}, 'red')
        map{i} = [(0:255)', zeros(256,2)]/255;
    elseif strcmp(mapstr{i}, 'green')
        map{i} = [zeros(256,1), (0:255)', zeros(256,1)]/255;
    elseif strcmp(mapstr{i}, 'blue')
        map{i} = [zeros(256,2), (0:255)']/255;
    elseif strcmp(mapstr{i}, 'gray')
        map{i} = gray(256);
    elseif strcmp(mapstr{i}, 'hot')
        map{i} = hot(256);
    elseif strcmp(mapstr{i}, 'hsv')
        map{i} = hsv(256);
    elseif strcmp(mapstr{i}, 'jet')
        map{i} = jet(256);
    elseif strcmp(mapstr{i}, 'parula')
        map{i} = parula(256);
    elseif strcmp(mapstr{i}, 'cool')
        map{i} = cool(256);
    elseif strcmp(mapstr{i}, 'spring')
        map{i} = spring(256);
    elseif strcmp(mapstr{i}, 'autumn')
        map{i} = autumn(256);
    elseif strcmp(mapstr{i}, 'summer')
        map{i} = summer(256);
    elseif strcmp(mapstr{i}, 'winter')
        map{i} = winter(256);
    elseif strcmp(mapstr{i}, 'bone')
        map{i} = bone(256);
    elseif strcmp(mapstr{i}, 'copper')
        map{i} = copper(256);
    elseif strcmp(mapstr{i}, 'pink')
        map{i} = pink(256);
    end
    I(i,:) = cellfun(@uint8, I(i,:), 'UniformOutput', 0);
    RGB0(i,:) = cellfun(@(x) ind2rgb(x, map{i}), I(i,:), 'UniformOutput', 0);
end
RGB = cellfun(@plus, RGB0(1,:), RGB0(2,:), 'UniformOutput', 0);
for i = 1 : size(RGB,2)
RGB{i}(RGB{i} > 1) = 1;
end
%handles.I = I;
%handles.map = map;
handles.RGB = RGB;
%handles.multiplier = mult;
guidata(hObject, handles);
ShowSlice(hObject, handles);

catch errorObj
    errordlg(getReport(errorObj,'extended','hyperlinks','off'),'Error');
end


function ShowSlice(~, handles)
%show the slice
p = handles.pixelsize;
%show = get(handles.popupmenu8, 'Value');
slice = str2double(get(handles.edit3,'String'));
% if show ~= 3
%     I = handles.I;
%     map = handles.map;
%     %imshow(I{show, slice}, map{show}, 'Parent', handles.axes1, 'XData', [0, size(I{show, slice},2) * p * 0.001], 'YData', [0, size(I{show, slice},1) * p * 0.001]);
%     image(ind2rgb(I{show, slice}, map{show}), 'Parent', handles.axes1);
%     set(handles.axes1,'XTickLabel',[], 'XTick',[], 'FontSize', 10);
% else
    RGB = handles.RGB;
    %imshow(RGB{slice}, 'Parent', handles.axes1, 'XData', [0, size(RGB{1,1},2) * p * 0.001], 'YData', [0, size(RGB{1,1},1) * p * 0.001]);
    image(RGB{slice}, 'Parent', handles.axes1, 'XData', [0, size(RGB{1,1},2) * p * 0.001], 'YData', [0, size(RGB{1,1},1) * p * 0.001]);
    set(handles.axes1,'XTickLabel',[], 'XTick',[], 'FontSize', 10);
% end


% --- Executes on selection change in popupmenu6.
function popupmenu6_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
updateimage(hObject, handles);



% --- Executes on selection change in popupmenu7.
function popupmenu7_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
updateimage(hObject, handles);




function edit6_Callback(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
updateimage(hObject, handles);




function edit7_Callback(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
updateimage(hObject, handles);




function edit8_Callback(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
GetSlices(hObject, handles);




function edit9_Callback(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
GetSlices(hObject, handles);




function GetSlices(hObject, handles)
% calculates Z-stack from eventlist
try
AB = handles.AB;
pixsizexy = str2double(get(handles.edit9, 'String'));
pixsizez = str2double(get(handles.edit8, 'String'));

fov = FOV(AB);
NZ = round(800/pixsizez); % number of slices in Z

I = cell(2,NZ); % I = {channel, Z}

for i = 1 : 2
if ~isempty (AB{i})
    I1 = slice( AB{i}, pixsizez, pixsizexy, fov );
    I(i,:) = I1;
else
    I(i,1:NZ) = {zeros(fov/pixsizexy)};
end
end
handles.I0 = I;
handles.pixelsize = pixsizexy;

% reset the current slice to 1
set(handles.edit3, 'String', num2str(1));

% set slider properties
slider_min = 1;
slider_max = NZ;
slider_step = [1/(slider_max-1), 1/(slider_max-1)];
slider_value = str2double(get(handles.edit3,'String'));
set(handles.slider1, 'Sliderstep', slider_step, 'Max', slider_max, 'Min', slider_min,'Value', slider_value);

guidata(hObject, handles);
updateimage(hObject, handles);

catch errorObj
    errordlg(getReport(errorObj,'extended','hyperlinks','off'),'Error');
end




function [Image] = slice (A, dz, pixsize, fov )
% cuts the eventlist on slices in Z and reconstructs an image in each
% Z-plane
% in nm
NZ = round(800/dz); % number of slices in Z
Image = cell (1, NZ);
l = size(A,1); % number of rows
n = ones(NZ,1);
Anew = cell(NZ,1);
for k = 1:l
    for i = 1:NZ
        b = i - NZ / 2;
    if (A(k,6) > dz * (b - 1)) && (A(k,6) <= dz * b)
        Anew{i}(n(i), :) = A(k,:);
        n(i) = n(i) + 1;
    end
    end
end
%parfor i = 1:NZ
for i = 1:NZ
    if ~isempty(Anew{i})
I = drpar(Anew{i}, pixsize, fov);
Image{i} = I;
    else
Image{i} = zeros(fov/pixsize);
    end
end


% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
