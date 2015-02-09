function [varargout] = contrastenh(varargin)
% CONTRASTENH MATLAB code for contrastenh.fig
%      CONTRASTENH, by itself, creates a new CONTRASTENH or raises the existing
%      singleton*.
%
%      H = CONTRASTENH returns the handle to a new CONTRASTENH or the handle to
%      the existing singleton*.
%
%      CONTRASTENH('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CONTRASTENH.M with the given input arguments.
%
%      CONTRASTENH('Property','Value',...) creates a new CONTRASTENH or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before contrastenh_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to contrastenh_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help contrastenh

% Last Modified by GUIDE v2.5 27-Jan-2015 17:52:11

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @contrastenh_OpeningFcn, ...
                   'gui_OutputFcn',  @contrastenh_OutputFcn, ...
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

% --- Executes just before contrastenh is made visible.
function contrastenh_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to contrastenh (see VARARGIN)

% Choose default command line output for contrastenh
handles.output = hObject;
handles.AB = varargin{1};
handles.ABcorr = varargin{2};
handles.folder = varargin{3};
handles.I1 = cell(2,3); % [redhist, redgauss, reddens; greenhist, greengauss, greendens]
handles.I0 = cell(2,3); % [redhist, redgauss, reddens; greenhist, greengauss, greendens]
handles.RGB1 = cell(3,1); % [hist, gauss, dens]
handles.RGB0 = cell(3,1); % [hist, gauss, dens]
handles.multiplier1 = cell(3,1); % [hist, gauss, dens]
handles.multiplier0 = cell(3,1); % [hist, gauss, dens]
handles.pixsize1 = zeros(3,1);
handles.pixsize0 = zeros(3,1);
handles.sigm = zeros(2,1);
% This sets up the initial plot - only do when we are invisible
% so window can get raised using contrastenh.
% Update handles structure
guidata(hObject, handles);
Reconstruct(hObject, handles);
% UIWAIT makes contrastenh wait for user response (see UIRESUME)
%uiwait(hObject);


% --- Outputs from this function are returned to the command line.
function varargout = contrastenh_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
varargout{2} = handles.RGB1;
varargout{3} = handles.multiplier1;
if isfield(handles, 'I0')
varargout{4} = handles.RGB0;
varargout{5} = handles.multiplier0;
end
%delete(hObject);


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
    % The GUI is no longer waiting, just close it
% end
RGB1 = handles.RGB1; % cell{hist, gauss}
RGB0 = handles.RGB0; % cell{hist, gauss}
multiplier1 = handles.multiplier1; % cell{hist, gauss}
multiplier0 = handles.multiplier0; % cell{hist, gauss}
saveorig = get(handles.checkbox3, 'Value'); %auto save original image
savecorr = get(handles.checkbox2, 'Value'); %auto save corrected image
pathname = handles.folder;

% auto save images in the current folder
if saveorig 
    if ~isempty(RGB0{1}) %hist
    filename = ['0h_' num2str(multiplier0{1}(1)) '_' num2str(multiplier0{1}(2)) '.tif'];
    FPName = [pathname filename];
    imwrite(RGB0{1}, FPName, 'Compression', 'lzw');
    end
    if ~isempty(RGB0{2}) %gauss
    filename = ['0g_' num2str(multiplier0{2}(1)) '_' num2str(multiplier0{2}(2)) '.tif'];
    FPName = [pathname filename];
    imwrite(RGB0{2}, FPName, 'Compression', 'lzw');
    end
    if ~isempty(RGB0{3}) %dens
    filename = ['0v_' num2str(multiplier0{3}(1)) '_' num2str(multiplier0{3}(2)) '.tif'];
    FPName = [pathname filename];
    imwrite(RGB0{3}, FPName, 'Compression', 'lzw');
    end
end
if savecorr
    if ~isempty(RGB1{1}) %hist
    filename = ['1h_' num2str(multiplier1{1}(1)) '_' num2str(multiplier1{1}(2)) '.tif'];
    FPName = [pathname filename];
    imwrite(RGB1{1}, FPName, 'Compression', 'lzw');
    end
    if ~isempty(RGB1{2}) %gauss
    filename = ['1g_' num2str(multiplier1{2}(1)) '_' num2str(multiplier1{2}(2)) '.tif'];
    FPName = [pathname filename];
    imwrite(RGB1{2}, FPName, 'Compression', 'lzw');
    end
    if ~isempty(RGB1{3}) %dens
    filename = ['1v_' num2str(multiplier1{3}(1)) '_' num2str(multiplier1{3}(2)) '.tif'];
    FPName = [pathname filename];
    imwrite(RGB1{3}, FPName, 'Compression', 'lzw');
    end
end

delete(hObject);



% --- Executes on selection change in popupmenu2.
function popupmenu2_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

UpdateImage(hObject, handles);



% --- Executes on selection change in popupmenu3.
function popupmenu3_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

UpdateImage(hObject, handles);




% --- Executes on selection change in popupmenu4.
function popupmenu4_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
UpdateImage(hObject, handles);


function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

UpdateImage(hObject, handles);



function edit4_Callback(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

UpdateImage(hObject, handles);



function UpdateImage(hObject, handles)

norm = get(handles.checkbox9, 'Value');
mult0(1) = str2double(get(handles.edit3, 'String'));
mult0(2) = str2double(get(handles.edit4, 'String'));
list = get(handles.radiobutton2, 'Value'); %1-corrected, 0-original;
mode = get(handles.popupmenu5, 'Value'); %1-histogram, 2-gaussian, 3-density;
if list
    I1 = handles.I1(:,mode);
elseif ~list
    I1 = handles.I0(:,mode);
end
cla(handles.axes1);
if ~isempty (I1{1}) || ~isempty (I1{2})
p = str2double(get(handles.edit5, 'String'));
show = get(handles.popupmenu4, 'Value');
value = get(handles.popupmenu2, 'Value'); 
stringlist = get(handles.popupmenu2, 'String');
mapstr{1} = stringlist{value}; %colormap for red eventlist
value = get(handles.popupmenu3, 'Value');
stringlist = get(handles.popupmenu3, 'String');
mapstr{2} = stringlist{value}; %colormap for green eventlist

mult = mult0;
if norm
    for k = 1:2
    if ~isempty (I1{k}) && max(max(I1{k})) > 0
        mult(k) = mult0(k) * 255 / max(max(I1{k}));
    end
    end
end

I = cell(2,1);
RGB0 = cell(2,1);
map = cell(2,1);
for i = 1 : 2
    I{i} = I1{i} * mult(i);
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
    I{i} = round(I{i});
    RGB0{i} = ind2rgb(I{i}, map{i});
end
RGB = RGB0{1} + RGB0{2};
RGB(RGB > 1) = 1;

iptsetpref('ImshowAxesVisible','on');
if show ~= 3
    imshow(I{show}, map{show}, 'Parent', handles.axes1, 'XData', [0, size(I{show},2) * p * 0.001], 'YData', [0, size(I{show},1) * p * 0.001]);
    set(handles.axes1,'XTickLabel',[], 'XTick',[], 'FontSize', 10);
else
   imshow(RGB, 'Parent', handles.axes1, 'XData', [0, size(I{1},2) * p * 0.001], 'YData', [0, size(I{1},1) * p * 0.001]);
   set(handles.axes1,'XTickLabel',[], 'XTick',[], 'FontSize', 10);
end
if list
    handles.RGB1{mode} = RGB;
    handles.multiplier1{mode} = mult;
elseif ~list
    handles.RGB0{mode} = RGB;
    handles.multiplier0{mode} = mult;
end
guidata(hObject, handles);

end





function edit5_Callback(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Reconstruct(hObject, handles);




% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Reconstruct(hObject, handles);


function Reconstruct(hObject, handles)
pixsize = str2double(get(handles.edit5, 'String'));
mode = get(handles.popupmenu5, 'Value');
% gaussmode = get(handles.radiobutton5, 'Value');
% histmode = get(handles.radiobutton4, 'Value');
draworig = get(handles.radiobutton3, 'Value');
drawcorr = get(handles.radiobutton2, 'Value');
sigma = str2double(get(handles.edit6, 'String'));
sigm = handles.sigm;

AB = handles.AB;
ABcorr = handles.ABcorr;
p0 = handles.pixsize0;
p1 = handles.pixsize1;

fov = FOV(AB);
zer = zeros(round(fov/pixsize));
Ic = cell(1,2);
I0 = cell(1,2);
if drawcorr
    for i = 1:2
        if ~isempty(ABcorr{i})
            if mode == 1 %hist
                if p1(1) ~= pixsize
                Ic{i} = draw(ABcorr{i}, pixsize, fov);
                end
            elseif mode == 2 %gauss
                if p1(2) ~= pixsize || sigm(1) ~= sigma
                    Ic{i} = gaussdraw(ABcorr{i}, pixsize, sigma, fov);
                end
            elseif mode == 3 %dens
                if p1(3) ~= pixsize
                    Ic{i} = drawVor(ABcorr{i}, pixsize, fov);
                end
            end
        else
            Ic{i} = zer;
        end
    end
elseif draworig
    for i = 1:2
        if ~isempty(AB{i})
            if mode == 1
                if p0(1) ~= pixsize
                I0{i} = draw(AB{i}, pixsize, fov);
                end
            elseif mode == 2
                if p0(2) ~= pixsize || sigm(2) ~= sigma
                    I0{i} = gaussdraw(AB{i}, pixsize, sigma, fov);
                end
            elseif mode == 3
                if p0(3) ~= pixsize
                    I0{i} = drawVor(AB{i}, pixsize, fov);
                end
            end
        else
            I0{i} = zer;
        end
    end
end
if mode == 2 && drawcorr
    if p1(2) ~= pixsize
        handles.I1(:,2) = Ic;
    end
        handles.pixsize1(2) = pixsize;
        handles.sigm(1) = sigma;
elseif mode == 1 && drawcorr
    if p1(1) ~= pixsize
        handles.I1(:,1) = Ic;
    end
        handles.pixsize1(1) = pixsize;
elseif mode == 3 && drawcorr
    if p1(3) ~= pixsize
        handles.I1(:,3) = Ic;
    end
        handles.pixsize1(3) = pixsize;
end
if mode == 2 && draworig
    if p0(2) ~= pixsize
        handles.I0(:,2) = I0;
    end
        handles.pixsize0(2) = pixsize;
        handles.sigm(2) = sigma;
elseif mode == 1 && draworig
    if p0(1) ~= pixsize
        handles.I0(:,1) = I0;
    end
        handles.pixsize0(1) = pixsize;
elseif mode == 3 && draworig
    if p0(3) ~= pixsize
        handles.I0(:,3) = I0;
    end
        handles.pixsize0(3) = pixsize;
end
guidata(hObject, handles);
UpdateImage(hObject, handles);


% --- Executes on button press in pushbutton7.
function pushbutton7_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
drawcorr = get(handles.radiobutton2, 'Value');
mode = get(handles.popupmenu5, 'Value');
if drawcorr
    RGB = handles.RGB1{mode};
elseif ~drawcorr
    RGB = handles.RGB0{mode};
end
if ~isempty(RGB)
pathname = handles.folder;
[FileName,PathName] = uiputfile({'*.tif'}, 'Save current image as', pathname);
if FileName ~= 0
FPName=[PathName FileName];
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

mode = get(handles.popupmenu5, 'Value');
% gaussmode = get(handles.radiobutton5, 'Value');
% histmode = get(handles.radiobutton4, 'Value');
draworig = get(handles.radiobutton3, 'Value');
drawcorr = get(handles.radiobutton2, 'Value');
p0 = handles.pixsize0;
p1 = handles.pixsize1;
if draworig && p0(mode) ~= 0
    set(handles.edit5, 'String', num2str(p0(mode)));
elseif drawcorr && p1(mode) ~= 0
    set(handles.edit5, 'String', num2str(p1(mode)));
end
guidata(hObject, handles);
Reconstruct(hObject, handles);


% --- Executes when selected object is changed in uipanel2.
function uipanel2_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel2 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
gaussmode = get(handles.radiobutton5, 'Value');
histmode = get(handles.radiobutton4, 'Value');
draworig = get(handles.radiobutton3, 'Value');
drawcorr = get(handles.radiobutton2, 'Value');
p0 = handles.pixsize0;
p1 = handles.pixsize1;
if gaussmode && draworig && p0(2) ~= 0
    set(handles.edit5, 'String', num2str(p0(2)));
elseif gaussmode && drawcorr && p1(2) ~= 0
    set(handles.edit5, 'String', num2str(p1(2)));
elseif histmode && draworig && p0(1) ~= 0
    set(handles.edit5, 'String', num2str(p0(1)));
elseif histmode && drawcorr && p1(1) ~= 0
    set(handles.edit5, 'String', num2str(p1(1)));
end
guidata(hObject, handles);
Reconstruct(hObject, handles);



function edit6_Callback(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Reconstruct(hObject, handles);


% --- Executes on selection change in popupmenu5.
function popupmenu5_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
mode = get(handles.popupmenu5, 'Value');
draworig = get(handles.radiobutton3, 'Value');
drawcorr = get(handles.radiobutton2, 'Value');
p0 = handles.pixsize0;
p1 = handles.pixsize1;
if mode == 2 && draworig && p0(2) ~= 0 %gauss
    set(handles.edit5, 'String', num2str(p0(2)));
elseif mode == 2 && drawcorr && p1(2) ~= 0 %gauss
    set(handles.edit5, 'String', num2str(p1(2)));
elseif mode == 1 && draworig && p0(1) ~= 0 %hist
    set(handles.edit5, 'String', num2str(p0(1)));
elseif mode == 1 && drawcorr && p1(1) ~= 0 %hist
    set(handles.edit5, 'String', num2str(p1(1)));
elseif mode == 3 && draworig && p0(3) ~= 0 %dens
    set(handles.edit5, 'String', num2str(p0(3)));
elseif mode == 3 && drawcorr && p1(3) ~= 0 %dens
    set(handles.edit5, 'String', num2str(p1(3)));
end
guidata(hObject, handles);
Reconstruct(hObject, handles);


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


% --- Executes on button press in checkbox9.
function checkbox9_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

UpdateImage(hObject, handles);
