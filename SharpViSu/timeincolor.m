function [varargout] = timeincolor(varargin)
% TIMEINCOLOR MATLAB code for timeincolor.fig
%      TIMEINCOLOR, by itself, creates a new TIMEINCOLOR or raises the existing
%      singleton*.
%
%      H = TIMEINCOLOR returns the handle to a new TIMEINCOLOR or the handle to
%      the existing singleton*.
%
%      TIMEINCOLOR('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TIMEINCOLOR.M with the given input arguments.
%
%      TIMEINCOLOR('Property','Value',...) creates a new TIMEINCOLOR or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before timeincolor_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to timeincolor_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help timeincolor

% Last Modified by GUIDE v2.5 23-Dec-2014 12:13:54

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @timeincolor_OpeningFcn, ...
                   'gui_OutputFcn',  @timeincolor_OutputFcn, ...
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

% --- Executes just before timeincolor is made visible.
function timeincolor_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to timeincolor (see VARARGIN)

% Choose default command line output for timeincolor
handles.output = hObject;
ABcorr = varargin{1};
ABorig = varargin{2};
handles.AB = cell(2,2); % [redcorr, redorig; greencorr, greenorig]
if size(ABcorr,1) == 1
    ABcorr = ABcorr.';
end
if size(ABorig,1) == 1
    ABorig = ABorig.';
end
handles.AB = [ABcorr, ABorig];
handles.folder = varargin{3};
handles.T = cell(2,2); % [redcorr, redorig; greencorr, greenorig]
handles.RGB = cell(2,2); % [redcorr, redorig; greencorr, greenorig]
handles.pixsize = zeros(2,2); % [redcorr, redorig; greencorr, greenorig]
handles.multiplier = zeros(2,2); % [redcorr, redorig; greencorr, greenorig]
% This sets up the initial plot - only do when we are invisible
% so window can get raised using timeincolor.
% Update handles structure
guidata(hObject, handles);
GetHSV(hObject, handles);


% --- Outputs from this function are returned to the command line.
function varargout = timeincolor_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
varargout{2} = handles.T;
varargout{3} = handles.RGB;


function GetHSV(hObject, handles)
try
if get(handles.radiobutton2, 'Value') && get(handles.radiobutton4, 'Value') % corr red
    index = [1, 1];
elseif get(handles.radiobutton2, 'Value') && ~get(handles.radiobutton4, 'Value') % corr green
    index = [2, 1];
elseif ~get(handles.radiobutton2, 'Value') && get(handles.radiobutton4, 'Value') % orig red
    index = [1, 2];
elseif ~get(handles.radiobutton2, 'Value') && ~get(handles.radiobutton4, 'Value') % orig green
    index = [2, 2];
end

A = handles.AB{index(1),index(2)};
if ~isempty(A)
p = str2double(get(handles.edit5, 'String'));
if p ~= handles.pixsize(index(1),index(2))

T0 = timeintHSV(A, p);
handles.T{index(1),index(2)} = T0;

handles.pixsize(index(1),index(2)) = p;

guidata(hObject, handles);
end
end
GetRGB(hObject, handles);

catch errorObj
    errordlg(getReport(errorObj,'extended','hyperlinks','off'),'Error');
end

function GetRGB(hObject, handles)
try
if get(handles.radiobutton2, 'Value') && get(handles.radiobutton4, 'Value') % corr red
    index = [1, 1];
elseif get(handles.radiobutton2, 'Value') && ~get(handles.radiobutton4, 'Value') % corr green
    index = [2, 1];
elseif ~get(handles.radiobutton2, 'Value') && get(handles.radiobutton4, 'Value') % orig red
    index = [1, 2];
elseif ~get(handles.radiobutton2, 'Value') && ~get(handles.radiobutton4, 'Value') % orig green
    index = [2, 2];
end

p = handles.pixsize(index(1),index(2));
multiplier = str2double(get(handles.edit4, 'String'));
if ~isempty(handles.T{index(1),index(2)})
    T = handles.T{index(1),index(2)};
    I = T;
    I(:,:,3) = multiplier * T(:,:,3);
    handles.RGB{index(1),index(2)} = hsv2rgb(I);
    handles.multiplier(index(1),index(2)) = multiplier;
end
cla(handles.axes1);
if ~isempty(handles.RGB{index(1),index(2)})
    RGB = handles.RGB{index(1),index(2)};
    iptsetpref('ImshowAxesVisible','on');
    imshow(RGB, 'parent', handles.axes1, 'XData', [0, size(RGB, 2) * p * 0.001], 'YData', [0, size(RGB, 1) * p * 0.001]);
    set(handles.axes1, 'XTickLabel', [], 'XTick', [], 'FontSize', 10);
end
guidata(hObject, handles);

catch errorObj
    errordlg(getReport(errorObj,'extended','hyperlinks','off'),'Error');
end



% --- Executes when user attempts to close figure.
function figure_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
% if isequal(get(hObject, 'waitstatus'), 'waiting')
%     % The GUI is still in UIWAIT, use UIRESUME
%     uiresume(hObject);
% else
%     The GUI is no longer waiting, just close it
% end
RGB = handles.RGB; 
saveorig = get(handles.checkbox3, 'Value'); %auto save original image
savecorr = get(handles.checkbox2, 'Value'); %auto save corrected image
pathname = handles.folder;

% auto save RGBs in the current folder
if saveorig 
    for i = 1:2
        if ~isempty(RGB{i,2})
            filename = 'T0_red.tif';
            FPName = [pathname filename];
            imwrite(RGB{i,2}, FPName, 'Compression', 'lzw');
        end
    end
end
if savecorr 
    for i = 1:2
        if ~isempty(RGB{i,1})
            filename = 'T1_red.tif';
            FPName = [pathname filename];
            imwrite(RGB{i,1}, FPName, 'Compression', 'lzw');
        end
    end
end

delete(hObject);





function edit4_Callback(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
GetRGB(hObject, handles);



function edit5_Callback(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
GetHSV(hObject, handles);



% --- Executes on button press in pushbutton7.
function pushbutton7_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if get(handles.radiobutton2, 'Value') && get(handles.radiobutton4, 'Value') % corr red
    index = [1, 1];
elseif get(handles.radiobutton2, 'Value') && ~get(handles.radiobutton4, 'Value') % corr green
    index = [2, 1];
elseif ~get(handles.radiobutton2, 'Value') && get(handles.radiobutton4, 'Value') % orig red
    index = [1, 2];
elseif ~get(handles.radiobutton2, 'Value') && ~get(handles.radiobutton4, 'Value') % orig green
    index = [2, 2];
end
if ~isempty(handles.RGB{index(1),index(2)})
    RGB = handles.RGB{index(1),index(2)};
    pathname = handles.folder;
    [FileName,PathName] = uiputfile({'*.tif'}, 'Save current image as', pathname);
    if FileName ~= 0
        FPName = [PathName FileName];
        imwrite(RGB, FPName, 'Compression', 'lzw');
    end
end


% --- Executes when selected object is changed in uipanel1.
function uipanel1_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel1 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
GetHSV(hObject, handles);



% --- Executes when selected object is changed in uipanel2.
function uipanel2_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel2 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
GetHSV(hObject, handles);
