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

% Last Modified by GUIDE v2.5 17-Mar-2015 11:09:57

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
handles.I1 = cell(2,5); % [redhist, redgauss, reddens; greenhist, greengauss, greendens]
handles.I0 = cell(2,5); % [redhist, redgauss, reddens; greenhist, greengauss, greendens]
handles.RGB1 = cell(5,1); % [hist, gauss, dens]
handles.RGB0 = cell(5,1); % [hist, gauss, dens]
handles.multiplier1 = cell(5,1); % [hist, gauss, dens]
handles.multiplier0 = cell(5,1); % [hist, gauss, dens]
handles.pixsize1 = zeros(5,1);
handles.pixsize0 = zeros(5,1);
handles.sigm = zeros(2,1);
handles.Voronoi = cell(2,5);
handles.Voronoi0 = cell(2,5);
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
    if ~isempty(RGB0{4}) %Qtree
    filename = ['0q_' num2str(multiplier0{4}(1)) '_' num2str(multiplier0{4}(2)) '.tif'];
    FPName = [pathname filename];
    imwrite(RGB0{4}, FPName, 'Compression', 'lzw');
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
    if ~isempty(RGB1{4}) %Q-Tree
    filename = ['1q_' num2str(multiplier1{4}(1)) '_' num2str(multiplier1{4}(2)) '.tif'];
    FPName = [pathname filename];
    imwrite(RGB1{4}, FPName, 'Compression', 'lzw');
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
try
monochr = get(handles.checkbox11, 'Value');
if monochr
    %set(handles.popupmenu2, 'Value', 3); % set string to gray
    set(handles.popupmenu2, 'Enable', 'off');
    %set(handles.popupmenu3, 'Value', 3); % set string to gray
    set(handles.popupmenu3, 'Enable', 'off');
else
    set(handles.popupmenu2, 'Enable', 'on');
    set(handles.popupmenu3, 'Enable', 'on');
end

norm = get(handles.checkbox9, 'Value');
mult0(1) = str2double(get(handles.edit3, 'String'));
mult0(2) = str2double(get(handles.edit4, 'String'));
list = get(handles.radiobutton2, 'Value'); %1-corrected, 0-original;
mode = get(handles.popupmenu5, 'Value'); %1-histogram, 2-gaussian, 3-density, 4-QuadTree;
if list
    I1 = handles.I1(:,mode);
    Voronoi = handles.Voronoi;
elseif ~list
    I1 = handles.I0(:,mode);
    Voronoi = handles.Voronoi0;
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
RGB0sh = cell(2,1);
map = cell(2,1);
for i = 1 : 2
    I{i} = I1{i} * mult(i);
    I{i} = round(I{i});
  if ~monochr
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
        RGB0{i} = ind2rgb(I{i}, map{i});
        RGB0sh{i} = RGB0{i};
  else
        map{i} = gray(256);
        RGB0{i} = I{i};
        RGB0sh{i} = ind2rgb(I{i}, map{i});
   end
end

if mode ~= 5
    
if ~isempty(RGB0{1}) && ~isempty(RGB0{2})
    RGB = RGB0{1} + RGB0{2};
    RGBsh = RGB0sh{1} + RGB0sh{2};
elseif ~isempty(RGB0{1})
    RGB = RGB0{1};
    RGBsh = RGB0sh{1};
elseif ~isempty(RGB0{2})
    RGB = RGB0{2};
    RGBsh = RGB0sh{2};
else
    RGB = zeros(900,900,3);
    RGBsh = zeros(900,900,3);
end
if ~monochr
    RGB(RGB > 1) = 1;
else
    RGB(RGB > 255) = 255;
end
RGBsh(RGBsh > 255) = 255;
iptsetpref('ImshowAxesVisible','on');

if show ~= 3
    imshow(I{show}, map{show}, 'Parent', handles.axes1, 'XData', [0, size(I{show},2) * p * 0.001], 'YData', [0, size(I{show},1) * p * 0.001]);
    set(handles.axes1,'XTickLabel',[], 'XTick',[], 'FontSize', 10);
else
   imshow(RGBsh, 'Parent', handles.axes1, 'XData', [0, size(I{1},2) * p * 0.001], 'YData', [0, size(I{1},1) * p * 0.001]);
   set(handles.axes1,'XTickLabel',[], 'XTick',[], 'FontSize', 10);
end
else            % show tessellations
    if show ~= 3
    RGB = RGB0{show};
    RGBsh = RGB0sh{show};
    

C = Voronoi{show,3};
V = Voronoi{show,2};
A = Voronoi{show,4};
hold on
    for i = 1:length(C)
%        if all(C{i} ~= 1)
            patch(V(C{i},1),V(C{i},2), RGBsh(i,:,:), 'Parent', handles.axes1);
%        end
    end
points = get(handles.checkbox12, 'Value');
if points
    plot(A(:,4),A(:,5),'.w', 'MarkerSize', 1);
end
hold off
xlim([0 18000]);
ylim([0 18000]);
set(gca, 'YDir', 'reverse');
    end
end
if list
    handles.RGB1{mode} = RGB;
    handles.multiplier1{mode} = mult;
elseif ~list
    handles.RGB0{mode} = RGB;
    handles.multiplier0{mode} = mult;
end
end
guidata(hObject, handles);

catch errorObj
    errordlg(getReport(errorObj,'extended','hyperlinks','off'),'Error');
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
try
pixsize = str2double(get(handles.edit5, 'String'));
mode = get(handles.popupmenu5, 'Value');
capacity = str2double(get(handles.edit7, 'String'));
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
Voronoi = handles.Voronoi;
Voronoi0 = handles.Voronoi0;

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
                    if isempty(Voronoi{i,1})
                         Voronoi(i,:) = VorArea(ABcorr{i});
                    end
                    Ic{i} = drawVor1(Voronoi{i,1}, Voronoi{i,4}, Voronoi{i,5}, pixsize, fov);
                end
            elseif mode == 4 %Quad-tree
                if p1(4) ~= capacity
                    Ic{i} = QuadTree(ABcorr{i}, capacity);
                end
            elseif mode == 5 %Voronoi decomposition
                if isempty(Voronoi{i,1})
                     Voronoi(i,:) = VorArea(ABcorr{i});
                end
                     Ic{i} = Voronoi{i,1};
            end
        else
            if any(mode == 1:3)
                Ic{i} = zer;
            else
                Ic{i} = [];
            end
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
                if isempty(Voronoi0{i,1})
                    Voronoi0(i,:) = VorArea(AB{i});
                end
                    I0{i} = drawVor1(Voronoi0{i,1}, Voronoi0{i,4}, Voronoi0{i,5}, pixsize, fov);
                end
            elseif mode == 4 %Quad-tree
                if p0(4) ~= capacity
                    cla(handles.axes1);
                    I0{i} = QuadTree(AB{i}, capacity);
                end
            elseif mode == 5 %Voronoi decomposition
                if isempty(Voronoi0{i,1})
                     Voronoi0(i,:) = VorArea(AB{i});
                end
                     I0{i} = Voronoi0{i,1};
            end
        else
            if any(mode == 1:3)
                I0{i} = zer;
            else
                I0{i} = [];
            end
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
elseif mode == 4 && drawcorr
    if p1(4) ~= capacity
        handles.I1(:,4) = Ic;
    end
        handles.pixsize1(4) = capacity;
elseif mode == 5 && drawcorr
        handles.I1(:,5) = Ic;
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
elseif mode == 4 && draworig
    if p0(4) ~= capacity
        handles.I0(:,4) = I0;
    end
        handles.pixsize0(4) = capacity;
elseif mode == 5 && draworig
%     if p0(4) ~= capacity
        handles.I0(:,5) = I0;
%     end
%         handles.pixsize0(4) = capacity;
end
handles.Voronoi = Voronoi;
handles.Voronoi0 = Voronoi0;
guidata(hObject, handles);
UpdateImage(hObject, handles);

catch errorObj
    errordlg(getReport(errorObj,'extended','hyperlinks','off'),'Error');
end


% --- Executes on button press in pushbutton7.
function pushbutton7_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try
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

catch errorObj
    errordlg(getReport(errorObj,'extended','hyperlinks','off'),'Error');
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
% if mode == 2
%     set(handles.edit6, 'Enable', 'on');
% else
%     set(handles.edit6, 'Enable', 'off');
% end
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
elseif mode == 4 && draworig && p0(4) ~= 0 %quadtree (capacity)
    set(handles.edit7, 'String', num2str(p0(4)));
elseif mode == 4 && drawcorr && p1(4) ~= 0 %quadtree (capacity)
    set(handles.edit7, 'String', num2str(p1(4)));
end
if mode == 5 %tessellations
    set(handles.checkbox12, 'Enable', 'on');
else
    set(handles.checkbox12, 'Enable', 'off');
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


% --- Executes on button press in checkbox11.
function checkbox11_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

UpdateImage(hObject, handles);



function edit7_Callback(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Reconstruct(hObject, handles);


% --- Executes during object creation, after setting all properties.
function edit7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu6.
function popupmenu6_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu6 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu6


% --- Executes during object creation, after setting all properties.
function popupmenu6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox12.
function checkbox12_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Reconstruct(hObject, handles);
