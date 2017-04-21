function varargout = aratea(varargin)
%LC Dec13 Incluyo el metodo de fast_oopsi para detectar las espigas. Dicho
%metodo calcula la probabilidad de espigas a partir de senhiales de calcio
%y depende de la tau y de la frecuencia de muestreo. 
%Me gusta mas el uso del gradiente ya que indica cuando habia valores en
%los transitorios de calcio que incrementaban con el tiempo, lo cual es de
%esperarse cuando la celula esta disparando.

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @aratea_OpeningFcn, ...
                   'gui_OutputFcn',  @aratea_OutputFcn, ...
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


% --- Executes just before aratea is made visible.
function aratea_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to aratea (see VARARGIN)

% Choose default command line output for aratea
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes aratea wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = aratea_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in newexp.
function newexp_Callback(hObject, eventdata, handles)
% hObject    handle to newexp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.video=uigetfile('*.*','Choose video');
coord_name=uigetfile('*.*','Choose coord');
coord_t=load(coord_name);
S=whos('coord_t');
if S.class=='struct'  %para cargar coordenadas en formato .mat o .txt LCsep15
coord_tf=struct2array(coord_t);
else
    coord_tf=coord_t;
end
handles.coord=coord_tf;

guidata(hObject, handles); %Guarda los handles nuevos


% --- Executes on button press in evalFFo.
function evalFFo_Callback(hObject, eventdata, handles)
% hObject    handle to evalFFo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.smooth=0; %Inicio de factores
handles.smoothfactor=0.15;
handles.threshold=50;
handles.cell_number=1;
set(handles.ncell, 'String', handles.cell_number);

a=FFo(handles.coord,handles.video);
b=FFoHalo(handles.coord,handles.video);
c=2*a-b;
handles.FFo_dat=c;

handles.cell_number=1;
set(handles.ncell, 'String', handles.cell_number);

%Grafica FFo para la celula #1; inicio de los parametros
axes(handles.axes1); 
cla;
plot(handles.FFo_dat(:,handles.cell_number))

%Grafica el gradiente y las espigas
axes(handles.axes2); 
cla;
dx=gradient(handles.FFo_dat(:,handles.cell_number));
% V.dt=1/2;P.gam=1-V.dt/5;dx=fast_oopsi(handles.FFo_dat(:,handles.cell_number),V,P);
std_dx=std(dx);
[size_dx,n_cells]=size(handles.FFo_dat);
hold
handles.del=ones(1,n_cells); %0 indica que la celula ha sido borrada
handles.spFFo_dat=zeros(size_dx,n_cells);
handles.spFFo_dat(:,handles.cell_number)=spikes(dx,handles.threshold,std_dx);
show(dx*handles.del(handles.cell_number),handles.threshold,std_dx,size_dx);
%borra las espigas de las celulas eliminadas
plot(dx)
hold
clc

guidata(hObject, handles); %Guarda los handles nuevos



% --- Executes on button press in cellant.
function cellant_Callback(hObject, eventdata, handles)
% hObject    handle to cellant (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

axes(handles.axes1);
cla;
if (handles.cell_number~=1)
handles.cell_number=handles.cell_number-1;
set(handles.ncell, 'String', handles.cell_number);
end

if (handles.smooth==0)
plot(handles.FFo_dat(:,handles.cell_number))
else
    plot(handles.smFFo_dat(:,handles.cell_number))
end


%Deriva y grafica abajo
axes(handles.axes2); 
cla;

if (handles.smooth==0)
dx=gradient(handles.FFo_dat(:,handles.cell_number)); 

% V.dt=1/2;P.gam=1-V.dt/5;dx=fast_oopsi(handles.FFo_dat(:,handles.cell_number),V,P);
std_dx=std(dx);
size_dx=size(handles.FFo_dat(:,handles.cell_number));
hold
handles.spFFo_dat(:,handles.cell_number)=spikes(dx,handles.threshold,std_dx);
show(dx*handles.del(handles.cell_number),handles.threshold,std_dx,size_dx);
else
dx=gradient(handles.smFFo_dat(:,handles.cell_number));
% V.dt=1/2;P.gam=1-V.dt/5;dx=fast_oopsi(handles.smFFo_dat(:,handles.cell_number),V,P);
std_dx=std(dx);
size_dx=size(handles.smFFo_dat(:,handles.cell_number));
hold
handles.spFFo_dat(:,handles.cell_number)=spikes(dx,handles.threshold,std_dx);
show(dx*handles.del(handles.cell_number),handles.threshold,std_dx,size_dx);
end

plot(dx)
hold
clc

guidata(hObject, handles); %Guarda los handles nuevos


% --- Executes on button press in cellnext.
function cellnext_Callback(hObject, eventdata, handles)
% hObject    handle to cellnext (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

axes(handles.axes1);
cla;
[size_dx,n_cells]=size(handles.FFo_dat);

if (handles.cell_number~=n_cells)
handles.cell_number=handles.cell_number+1;
set(handles.ncell, 'String', handles.cell_number);
end

if (handles.smooth==0)
plot(handles.FFo_dat(:,handles.cell_number))
else
    plot(handles.smFFo_dat(:,handles.cell_number))
end


%Deriva y grafica abajo
axes(handles.axes2); 
cla;

if (handles.smooth==0)
dx=gradient(handles.FFo_dat(:,handles.cell_number));
% V.dt=1/2;P.gam=1-V.dt/5;dx=fast_oopsi(handles.FFo_dat(:,handles.cell_number),V,P);

std_dx=std(dx);
size_dx=size(handles.FFo_dat(:,handles.cell_number));
hold
handles.spFFo_dat(:,handles.cell_number)=spikes(dx,handles.threshold,std_dx);
show(dx*handles.del(handles.cell_number),handles.threshold,std_dx,size_dx);
else
dx=gradient(handles.smFFo_dat(:,handles.cell_number));
% V.dt=1/2;P.gam=1-V.dt/5;dx=fast_oopsi(handles.smFFo_dat(:,handles.cell_number),V,P);
std_dx=std(dx);
size_dx=size(handles.smFFo_dat(:,handles.cell_number));
hold
handles.spFFo_dat(:,handles.cell_number)=spikes(dx,handles.threshold,std_dx);
show(dx*handles.del(handles.cell_number),handles.threshold,std_dx,size_dx);
end

plot(dx)
hold
clc

guidata(hObject, handles); %Guarda los handles nuevos



function smoothfactor_Callback(hObject, eventdata, handles)
% hObject    handle to smoothfactor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of smoothfactor as text
span=str2double(get(hObject,'String')); % returns contents of smoothfactor as a double

if isnan(span)
    set(hObject, 'String', 0.15);
    errordlg('Input must be a number','Error');
end

handles.smoothfactor=span;

% handles.span=span;
guidata(hObject, handles); %Guarda los handles nuevos



% --- Executes during object creation, after setting all properties.
function smoothfactor_CreateFcn(hObject, eventdata, handles)
% hObject    handle to smoothfactor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on button press in smooth/filter.
function smooth_Callback(hObject, eventdata, handles)
% hObject    handle to smooth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

axes(handles.axes1);
cla;
if (handles.smooth==0)
    handles.smooth=handles.smooth + 1;
end

% if (handles.smoothfactor==0)
%     handles.span=2;
%     handles.smooth=handles.smooth + 1;
%     handles.smoothfactor=1;
%     span=handles.span;
% elseif (handles.smooth==0)
%     span = handles.span; % Size of the averaging window
%     handles.smooth=handles.smooth + 1;
% else
%     span = handles.span*handles.smooth; % Size of the averaging window
%     handles.smooth=handles.smooth + 1;
% end

%window = ones(span,1)/span;
%handles.ssFFo_dat(:,handles.cell_number)=convn(handles.sFFo_dat(:,handles.cell_number),window,'valid');

% handles.smFFo_dat=convn(handles.FFo_dat,window,'valid'); %'same' en vez de 'valid' para que el arreglo de salida sea igual
%LCDec13 cambio de filtro para no disminuir el tamanhio del arreglo uso un filtro digital en
%vez de la convolucion en n
%handles.smFFo_dat=filter(ones(1,span)/span,1,handles.FFo_dat); %LC Dec13
%Uso el filto paso bajas de ct2
handles.smFFo_dat=lowpass_filter_ct2(double(handles.FFo_dat'),handles.smoothfactor)';
%LCEn14 mantiene el mismo numero de puntos uso una ventana del 10% total de
%puntos considerando una frecuencia de muestreo ~4 Hz. Funciona relindo


plot(handles.smFFo_dat(:,handles.cell_number))

%Derivada
handles.dFFo_dat=gradient(handles.smFFo_dat);
% V.dt=1/2;P.gam=1-V.dt/5;handles.dFFo_dat=fast_oopsi(handles.smFFo_dat(:,handles.cell_number),V,P);

%Deriva y grafica abajo
axes(handles.axes2); 
cla;
dx=gradient(handles.smFFo_dat(:,handles.cell_number));
% V.dt=1/2;P.gam=1-V.dt/5;dx=fast_oopsi(handles.smFFo_dat(:,handles.cell_number),V,P);

std_dx=std(dx);
[size_dx,n_cells]=size(handles.smFFo_dat);
hold
handles.spFFo_dat=zeros(size_dx,n_cells);
handles.spFFo_dat(:,handles.cell_number)=spikes(dx,handles.threshold,std_dx);
show(dx*handles.del(handles.cell_number),handles.threshold,std_dx,size_dx);
plot(dx)
hold
clc

guidata(hObject, handles); %Guarda los handles nuevos



% --- Executes on button press in reset.
function reset_Callback(hObject, eventdata, handles)
% hObject    handle to reset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.smooth=0;
% handles.smoothfactor=0.15;
% handles.threshold=50;

axes(handles.axes1);
cla;

plot(handles.FFo_dat(:,handles.cell_number))

%Deriva y grafica abajo
axes(handles.axes2); 
cla;

[size_dx,n_cells]=size(handles.FFo_dat);
hold
handles.spFFo_dat=zeros(size_dx,n_cells);
handles.dFFo_dat=zeros(size_dx,n_cells);

for ii=1:n_cells
   handles.dFFo_dat(:,ii)=gradient(handles.FFo_dat(:,ii));
%    V.dt=1/2;P.gam=1-V.dt/5; handles.dFFo_dat(:,ii)=fast_oopsi(handles.FFo_dat(:,ii),V,P);

    dxsp=handles.dFFo_dat(:,ii);
    std_dx=std(dxsp);
    handles.spFFo_dat(:,ii)=spikes(dxsp,handles.threshold,std_dx);    
end

dx=gradient(handles.FFo_dat(:,handles.cell_number));
% V.dt=1/2;P.gam=1-V.dt/5;dx=fast_oopsi(handles.FFo_dat(:,handles.cell_number),V,P);

std_dx=std(dx);
show(dx*handles.del(handles.cell_number),handles.threshold,std_dx,size_dx);
plot(dx)
hold
clc

handles.del=ones(1,n_cells); %0 indica que la celula ha sido borrada

guidata(hObject, handles); %Guarda los handles nuevos





% --- Executes on button press in savedata.
function savedata_Callback(hObject, eventdata, handles)
% hObject    handle to savedata (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Derivada y espigas de smFFo_dat 
if (handles.smooth==0)
   handles.smFFo_dat=handles.FFo_dat;
end
    
    
[size_dx,n_cells]=size(handles.smFFo_dat);
handles.dFFo_dat=zeros(size_dx,n_cells);
handles.spFFo_dat=zeros(size_dx,n_cells);


%Renumera las celulas activas y elimina las borradas y las no activas
jj=1;
for ii=1:n_cells
    if (handles.del(ii)==1) 
    handles.dFFo_dat(:,ii)=gradient(handles.smFFo_dat(:,ii));
%    V.dt=1/2;P.gam=1-V.dt/5; handles.dFFo_dat(:,ii)=fast_oopsi(handles.smFFo_dat(:,ii),V,P);

    dx=handles.dFFo_dat(:,ii);
    std_dx=std(dx);
    handles.spFFo_dat(:,ii)=spikes(dx,handles.threshold,std_dx);
       if (sum(handles.spFFo_dat(:,ii))>=1)
       handles.renFFo_dat(:,jj)=handles.FFo_dat(:,ii);    
       handles.rensmFFo_dat(:,jj)=handles.smFFo_dat(:,ii);
       handles.rendFFo_dat(:,jj)=handles.dFFo_dat(:,ii);
       handles.renspFFo_dat(:,jj)=handles.spFFo_dat(:,ii);
       handles.rencoord(jj,:)=handles.coord(ii,:);
       jj=jj+1;
       end
    end
end
    
handles.renspplotFFo_dat=RasNum(handles.renspFFo_dat);

figure
plot(handles.rensmFFo_dat)

% figure
% plot(handles.renspplotFFo_dat,'+')

%Salva los archivos como .dat
FFo=handles.renFFo_dat;
save FFo.dat FFo -ascii
save FFo FFo

FFo_smooth=handles.rensmFFo_dat;
save FFo_smooth.dat FFo_smooth -ascii
save FFo_smooth FFo_smooth 

FFo_deriv=handles.rendFFo_dat;
save FFo_deriv.dat FFo_deriv -ascii
save FFo_deriv FFo_deriv

Z=handles.renspFFo_dat;
Spikes=Z';                     %Para que quede una matriz [cells,frames]
save Spikes.dat Spikes -ascii
save Spikes Spikes

% W=handles.renspplotFFo_dat;
% save Spikes_plot.dat W -ascii

Coord_active=handles.rencoord;
save Coord_active.dat Coord_active -ascii
save Coord_active Coord_active

% U=trail(Zb);
% save Spikes_trail.dat U -ascii
% 
% T=RasNum(U');
% T=T';
% save Spikes_trail_plot.dat T -ascii

figure
% plot(T','+')
imagesc(Spikes);

guidata(hObject, handles); %Guarda los handles nuevos




% --- Executes on button press in deletecell.
function deletecell_Callback(hObject, eventdata, handles)
% hObject    handle to deletecell (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.spFFo_dat(:,handles.cell_number)=zeros();
handles.del(handles.cell_number)=0; %Borra la celula seleccionada

%Deriva y grafica abajo
axes(handles.axes2); 
cla;

guidata(hObject, handles); %Guarda los handles nuevos

function threshold_Callback(hObject, eventdata, handles)
% hObject    handle to threshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of threshold as text
%        str2double(get(hObject,'String')) returns contents of threshold as a double

umbral=str2double(get(hObject,'String')); % returns contents of smoothfactor as a double

if isnan(umbral)
    set(hObject, 'String', 50);
    errordlg('Input must be a number','Error');
end

handles.threshold=umbral;

% handles.span=span;
guidata(hObject, handles); %Guarda los handles nuevos

% --- Executes during object creation, after setting all properties.
function threshold_CreateFcn(hObject, eventdata, handles)
% hObject    handle to threshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
