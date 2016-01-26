function varargout = Zcalibration(varargin)
% ZCALIBRATION MATLAB code for Zcalibration.fig
%      ZCALIBRATION, by itself, creates a new ZCALIBRATION or raises the existing
%      singleton*.
%
%      H = ZCALIBRATION returns the handle to a new ZCALIBRATION or the handle to
%      the existing singleton*.
%
%      ZCALIBRATION('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ZCALIBRATION.M with the given input arguments.
%
%      ZCALIBRATION('Property','Value',...) creates a new ZCALIBRATION or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Zcalibration_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Zcalibration_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Zcalibration

% Last Modified by GUIDE v2.5 21-Nov-2014 14:32:52

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Zcalibration_OpeningFcn, ...
                   'gui_OutputFcn',  @Zcalibration_OutputFcn, ...
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


% --- Executes just before Zcalibration is made visible.
function Zcalibration_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Zcalibration (see VARARGIN)

% Choose default command line output for Zcalibration
handles.output = hObject;
handles.ABC = cell(3,1);
if nargin == 1
handles.folder = varargin{1};
end
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Zcalibration wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Zcalibration_OutputFcn(hObject, eventdata, handles) 
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


function GetFile( hObject, handles, channel)
try
format = get(handles.popupmenu1, 'Value');
    extension = '*.*';
if isfield(handles, 'folder')
    pathname = handles.folder;
    [filename, pathname] = uigetfile({extension}, 'Please select eventlists of the Z-stack', pathname, 'MultiSelect', 'on');
else
    [filename, pathname] = uigetfile({extension}, 'MultiSelect', 'on');
end
if pathname ~= 0
    if channel == 1 % 642
        set (handles.edit2, 'String', pathname);
    elseif channel == 2 % 532
        set (handles.edit4, 'String', pathname);
    elseif channel == 3 % 488
        set (handles.edit5, 'String', pathname);
    end
if ~iscell(filename)
    filenamet = cell(1);
    filenamet{1} = filename;
    filename = filenamet;
end
A = cell(size(filename, 2),1);
for i = 1:size(filename, 2)
    B = importdata([pathname filename{i}]);
    if isstruct(B)
        B = B.data;
    end
A{i} = adjustformat(B, format);    % adjust the format to the standard
end
handles.folder = pathname;
handles.ABC{channel} = A;
else
    if channel == 1 % 642
        set (handles.edit2, 'String', []);
    elseif channel == 2 % 532
        set (handles.edit4, 'String', []);
    elseif channel == 3 % 488
        set (handles.edit5, 'String', []);
    end
    handles.ABC{channel} = [];
end

cla(handles.axes1); %clear the fit plot
guidata(hObject,handles);

catch errorObj
    errordlg(getReport(errorObj,'extended','hyperlinks','off'),'Error');
end




% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
GetFile(hObject, handles, 2);



% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
GetFile(hObject, handles, 3);



% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try
ABC = handles.ABC;
step = str2double(get(handles.edit6, 'String'));
degree = str2double(get(handles.edit10, 'String'));
absZ = get(handles.checkbox1, 'Value');
p = cell(3,1);
A = cell(3,1);
An = cell(3,1);
Am = cell(3,1);
Z0 = zeros(3,1);
STD = cell(3,1);
for i = 1:3
if ~isempty(ABC{i})
B = ABC{i};
A{i} = zeros(size(B,1),4);
An{i} = A{i};
if i == 1
Zoffset{i} = str2double(get(handles.edit7, 'String'));
Soffset{i} = str2double(get(handles.edit14, 'String'));
elseif i == 2
Zoffset{i} = str2double(get(handles.edit8, 'String'));
Soffset{i} = str2double(get(handles.edit15, 'String'));
elseif i == 3
Zoffset{i} = str2double(get(handles.edit9, 'String'));
Soffset{i} = str2double(get(handles.edit16, 'String'));
end
STD{i} = zeros(size(B,1),1);
for j = 1:size(B,1)
    C = B{j};
    At = zeros(size(C,1),2);
    At(:,1) = step * (j-1);
    At(:,2) = C(:,8) - C(:,9);
    STD{i}(j) = std(At(:,2));
    Am{i} = [Am{i}; At];
mX = mean(C(:, 8));
mY = mean(C(:, 9));
A{i}(j, 1) = step * (j-1); % Z position
A{i}(j, 2) = mX; % mean sigmaX
A{i}(j, 3) = mY; % mean sigmaY
A{i}(j, 4) = mX - mY; % difference 
end
p0 = polyfit(A{i}(:,4), A{i}(:,1), degree);
if get(handles.checkbox2, 'Value') % if center Z
Z0(i) = polyval(p0,0); % evaluate Z(dsimga) at sigma = 0
end
if ~absZ
An{i} = A{i};
An{i}(:,4) = A{i}(:,4) - Soffset{i};
An{i}(:,1) = A{i}(:,1) - Z0(i) - Zoffset{i};

p{i} = polyfit(An{i}(:,4), An{i}(:,1), degree);

Am{i}(:,2) = Am{i}(:,2) - Soffset{i};
Am{i}(:,1) = Am{i}(:,1) - Z0(i) - Zoffset{i};
end
end
end

if absZ
    Z0(Z0 == 0) = [];
    Zconst = mean(Z0);
    if isnan(Zconst)
        Zconst = 0;
    end
    for i = 1:3
        if ~isempty(ABC{i})
            An{i} = A{i};
            An{i}(:,4) = A{i}(:,4);
            An{i}(:,1) = A{i}(:,1) - Zconst;

            p{i} = polyfit(An{i}(:,4), An{i}(:,1), degree);
                        
            Am{i}(:,1) = Am{i}(:,1) - Zconst;
        end
    end
end
handles.STD = STD;
handles.p = p;
handles.A = A;
handles.An = An;
handles.Am = Am;
guidata(hObject,handles);
UpdatePlot(hObject, handles);

if isempty(ABC{1}) && isempty(ABC{2}) && isempty(ABC{3})
    h = errordlg('Please load calibration data first','No calibration data');
end
catch errorObj
    errordlg(getReport(errorObj,'extended','hyperlinks','off'),'Error');
end

function UpdatePlot(~, handles)
cla(handles.axes1); %clear the fit plot
An = handles.An;
Am = handles.Am;
STD = handles.STD;
p = handles.p;
ShowAll = get(handles.checkbox3, 'Value');
graph = get(handles.popupmenu2, 'Value');
color{1} = 'r'; color{2} = 'g'; color{3} = 'b'; 
if ShowAll
    
if graph == 1 % Z(deltaSigma)
    set(handles.text11, 'String', 'sx - sy, nm');
    set(handles.text12, 'String', 'Z, nm');
hold on
for i = 1:3
    if ~isempty(An{i})
plot(handles.axes1, Am{i}(:,2), Am{i}(:,1), ['.' color{i}], An{i}(:,4), polyval(p{i}, An{i}(:,4)), ['-' color{i}], 'MarkerSize', 6); % X = sigmaX-sigmaY, nm; Y = Z, nm
    end
end
hold off
elseif graph == 2 % deltaSigma(Z)
    set(handles.text11, 'String', 'Z, nm');
    set(handles.text12, 'String', 'sx-sy, nm');
hold on
for i = 1:3
    if ~isempty(An{i})
        plot(handles.axes1, Am{i}(:,1), Am{i}(:,2), ['.' color{i}], polyval(p{i}, An{i}(:,4)), An{i}(:,4), ['-' color{i}], 'MarkerSize', 6); % X = Z, nm; Y = sigmaX-sigmaY, nm
        
    end
end
hold off    
end

else % show only average experimental points
    
if graph == 1 % Z(deltaSigma)
    set(handles.text11, 'String', 'sx - sy, nm');
    set(handles.text12, 'String', 'Z, nm');
hold on
for i = 1:3
    if ~isempty(An{i})
plot(handles.axes1, An{i}(:,4), An{i}(:,1), ['.' color{i}], An{i}(:,4), polyval(p{i}, An{i}(:,4)), ['-' color{i}]); % X = sigmaX-sigmaY, nm; Y = Z, nm

    end
end
hold off
elseif graph == 2 % deltaSigma(Z)
    set(handles.text11, 'String', 'Z, nm');
    set(handles.text12, 'String', 'sx-sy, nm');
hold on
for i = 1:3
    if ~isempty(An{i})
plot(handles.axes1, An{i}(:,1), An{i}(:,4), ['.' color{i}], polyval(p{i}, An{i}(:,4)), An{i}(:,4), ['-' color{i}]); % X = Z, nm; Y = sigmaX-sigmaY, nm
errorbar(An{i}(:,1), An{i}(:,4), STD{i}, ['.' color{i}], 'MarkerSize', 6, 'Parent', handles.axes1);
    end
end
hold off    
end

end


% % save the plot
% 
% pathname = handles.folder;
% FPName = [pathname 'Calibration_curve.png'];
% h2 = figure('visible', 'off');
% ax = axes;
% 
% hold on
% for i = 1:3
%     if ~isempty(An{i})
% plot(ax, An{i}(:,1), An{i}(:,4), ['.' color{i}], polyval(p{i}, -330:10:180), -330:10:180, ['-' color{i}], 'LineWidth', 1.5, 'MarkerSize', 12); % X = Z, nm; Y = sigmaX-sigmaY, nm
% errorbar(An{i}(:,1), An{i}(:,4), STD{i}, ['.' color{i}], 'MarkerSize', 4, 'LineWidth', 1.5, 'Parent', ax);
%     end
% end
% hold off    
% axis([-300, 400, -350, 200]);
% xlabel('Z, nm', 'FontSize', 18);
% ylabel('sigmax - sigmay, nm', 'FontSize', 18);
% set(gca, 'FontSize', 18, 'LineWidth', 2);
% myStyle = hgexport('factorystyle');
% myStyle.Format = 'png';
% myStyle.Width = 9;
% myStyle.Height = 6;
% myStyle.Resolution = 300;
% myStyle.Units = 'inch';
% myStyle.FixedFontSize = 8;
% hgexport(h2, FPName, myStyle, 'Format', 'png');
% close(h2);



% --- Executes on button press in pushbutton7.
function pushbutton7_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
p = handles.p;
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
FileName{1} = '\Z642.dat';
FileName{2} = '\Z532.dat';
FileName{3} = '\Z488.dat';
for i = 1 : 3
    if ~isempty(p{i})
FPName=[currentDir FileName{i}];
dlmwrite(FPName, p{i});
    end
end

if isempty(p{1}) && isempty(p{2}) && isempty(p{3})
    h = errordlg('Please fit the data first','No fit');
end


% --- Executes on button press in checkbox1.
function checkbox1_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
absZ = get(hObject,'Value');
if absZ
    set(handles.edit7, 'String', 0, 'Style', 'text', 'BackgroundColor', [0.94 0.94 0.94], 'ForegroundColor', [0.5 0.5 0.5]);
    set(handles.edit8, 'String', 0, 'Style', 'text', 'BackgroundColor', [0.94 0.94 0.94], 'ForegroundColor', [0.5 0.5 0.5]);
    set(handles.edit9, 'String', 0, 'Style', 'text', 'BackgroundColor', [0.94 0.94 0.94], 'ForegroundColor', [0.5 0.5 0.5]);
    set(handles.edit14, 'String', 0, 'Style', 'text', 'BackgroundColor', [0.94 0.94 0.94], 'ForegroundColor', [0.5 0.5 0.5]);
    set(handles.edit15, 'String', 0, 'Style', 'text', 'BackgroundColor', [0.94 0.94 0.94], 'ForegroundColor', [0.5 0.5 0.5]);
    set(handles.edit16, 'String', 0, 'Style', 'text', 'BackgroundColor', [0.94 0.94 0.94], 'ForegroundColor', [0.5 0.5 0.5]);
else
    set(handles.edit7, 'Style', 'edit', 'BackgroundColor', [1 1 1], 'ForegroundColor', [0 0 0]);
    set(handles.edit8, 'Style', 'edit', 'BackgroundColor', [1 1 1], 'ForegroundColor', [0 0 0]);
    set(handles.edit9, 'Style', 'edit', 'BackgroundColor', [1 1 1], 'ForegroundColor', [0 0 0]);
    set(handles.edit14, 'Style', 'edit', 'BackgroundColor', [1 1 1], 'ForegroundColor', [0 0 0]);
    set(handles.edit15, 'Style', 'edit', 'BackgroundColor', [1 1 1], 'ForegroundColor', [0 0 0]);
    set(handles.edit16, 'Style', 'edit', 'BackgroundColor', [1 1 1], 'ForegroundColor', [0 0 0]);
end



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
