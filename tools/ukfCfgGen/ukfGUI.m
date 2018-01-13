function varargout = ukfGUI(varargin)
% UKFGUI MATLAB code for ukfGUI.fig
%      UKFGUI, by itself, creates a new UKFGUI or raises the existing
%      singleton*.
%
%      H = UKFGUI returns the handle to a new UKFGUI or the handle to
%      the existing singleton*.
%
%      UKFGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in UKFGUI.M with the given input
%      arguments.
%
%      UKFGUI('Property','Value',...) creates a new UKFGUI or raises
%      the existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ukfGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ukfGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ukfGUI

% Last Modified by GUIDE v2.5 13-Jan-2018 22:53:49

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ukfGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @ukfGUI_OutputFcn, ...
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

% --- Executes just before ukfGUI is made visible.
function ukfGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ukfGUI (see VARARGIN)

% Choose default command line output for ukfGUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

initialize_gui(hObject, handles, false);

% UIWAIT makes ukfGUI wait for user response (see UIRESUME)
% uiwait(handles.UKF_CFG_GEN);


% --- Outputs from this function are returned to the command line.
function varargout = ukfGUI_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
handles.output = hObject;
handles = guidata(hObject);
varargout{1} = handles.output;


% --- Executes during object creation, after setting all properties.
function Pxx_value_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Pxx_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Pxx_value_Callback(hObject, eventdata, handles)
% hObject    handle to Pxx_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Pxx_value as text
%        str2double(get(hObject,'String')) returns contents of Pxx_value as a double
Pxx = eval(get(hObject, 'String'));
if isnan(Pxx)
    set(hObject, 'String', 0);
    errordlg('Input must be a number','Error');
end

% Save the new Pxx_value value
handles.ukfdata.Pxx = Pxx;
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function Ryy_value_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Ryy_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Ryy_value_Callback(hObject, eventdata, handles)
% hObject    handle to Ryy_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Ryy_value as text
%        str2double(get(hObject,'String')) returns contents of Ryy_value as a double
Ryy = eval(get(hObject, 'String'));
if isnan(Ryy)
    set(hObject, 'String', 0);
    errordlg('Input must be a number','Error');
end

% Save the new Ryy_value value
handles.ukfdata.Ryy = Ryy;
guidata(hObject,handles)

% --- Executes on button press in Generate_ukfCfg.
function Generate_ukfCfg_Callback(hObject, eventdata, handles)
% hObject    handle to Generate_ukfCfg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ukfCfgGen(handles);

% --- Executes on button press in reset.
function reset_Callback(hObject, eventdata, handles)
% hObject    handle to reset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

initialize_gui(gcbf, handles, true);

% --- Executes when selected object changed in unitgroup.
function unitgroup_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in unitgroup 
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if (hObject == handles.english)
    set(handles.text4, 'String', 'lb/cu.in');
    set(handles.text5, 'String', 'cu.in');
    set(handles.text6, 'String', 'lb');
else
    set(handles.text4, 'String', 'kg/cu.m');
    set(handles.text5, 'String', 'cu.m');
    set(handles.text6, 'String', 'kg');
end

% --------------------------------------------------------------------
function initialize_gui(fig_handle, handles, isreset)
% If the ukfdata field is present and the reset flag is false, it means
% we are we are just re-initializing a GUI by calling it from the cmd line
% while it is up. So, bail out as we dont want to reset the data.
if isfield(handles, 'ukfdata') && ~isreset
    return;
end
handles.ukfdata.StateFcn = '{''x(1) = x(1) + dT*x(2);'';''x(2) = (1 - dT*0.1)*x(2) - dT*16.003263*sin(x(1));''}';
set(handles.StateFcn,  'String', handles.ukfdata.StateFcn);

handles.ukfdata.MeasFcn = '{''y(1) = x(1);''}';
set(handles.MeasFcn,  'String', handles.ukfdata.MeasFcn);

[xL,~] = size(handles.ukfdata.StateFcn);
[yL,~] = size(handles.ukfdata.MeasFcn);
sL = 2*xL+1;
uL= 2;

handles.ukfdata.Pxx = diag(ones(1,xL)*0.01);
handles.ukfdata.Qxx = diag(ones(1,xL)*0.2);
handles.ukfdata.Ryy  = 0.13;
handles.ukfdata.dT  = 0.0001;
handles.ukfdata.cfgid  = 1;

handles.ukfdata.alpha = 0.1;
set(handles.alphaEdit,  'String', handles.ukfdata.alpha);
set(handles.alphaSlider,  'Value', handles.ukfdata.alpha);

handles.ukfdata.betha = 2;
set(handles.bethaEdit,  'String', handles.ukfdata.betha);
set(handles.bethaSlider,  'Value', handles.ukfdata.betha);

handles.ukfdata.kappa = 0;
set(handles.kappaEdit,  'String', handles.ukfdata.kappa);
set(handles.kappaSlider,  'Value', handles.ukfdata.kappa);

handles.ukfdata.LimitsEnable = false;

%set(handles.Pxx_value, 'String', handles.ukfdata.Pxx);
%set(handles.Ryy_value,  'String', handles.ukfdata.Ryy);
%set(handles.Qxx_value,  'String', handles.ukfdata.Ryy);

% Update handles structure
%guidata(handles.UKF_CFG_GEN, handles);
guidata(fig_handle,handles)



function Qxx_value_Callback(hObject, eventdata, handles)
% hObject    handle to Qxx_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Qxx_value as text
%        str2double(get(hObject,'String')) returns contents of Qxx_value as a double
Qxx = eval(get(hObject, 'String'));
if isnan(Qxx)
    set(hObject, 'String', 0);
    errordlg('Input must be a number','Error');
end

% Save the new Ryy_value value
handles.ukfdata.Qxx = Qxx;
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function Qxx_value_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Qxx_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit7_Callback(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit7 as text
%        str2double(get(hObject,'String')) returns contents of edit7 as a double


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


% --- Executes on button press in LimitsEnable.
function LimitsEnable_Callback(hObject, eventdata, handles)
% hObject    handle to LimitsEnable (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of LimitsEnable
handles.ukfdata.LimitsEnable = get(hObject,'Value');
guidata(hObject,handles);


function StateTransitionFcn_Callback(hObject, eventdata, handles)
% hObject    handle to StateTransitionFcn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of StateTransitionFcn as text
%        str2double(get(hObject,'String')) returns contents of StateTransitionFcn as a double
% StateFcn = eval(get(hObject, 'String'));
% if isnan(StateFcn)
%     set(hObject, 'String', 0);
%     errordlg('Input must be a number','Error');
% end
% 
% % Save the new Ryy_value value
% handles.ukfdata.StateFcn = StateFcn;
% guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function StateTransitionFcn_CreateFcn(hObject, eventdata, handles)
% hObject    handle to StateTransitionFcn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit9_Callback(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit9 as text
%        str2double(get(hObject,'String')) returns contents of edit9 as a double


% --- Executes during object creation, after setting all properties.
function edit9_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit10_Callback(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit10 as text
%        str2double(get(hObject,'String')) returns contents of edit10 as a double


% --- Executes during object creation, after setting all properties.
function edit10_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit11_Callback(hObject, eventdata, handles)
% hObject    handle to edit11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit11 as text
%        str2double(get(hObject,'String')) returns contents of edit11 as a double


% --- Executes during object creation, after setting all properties.
function edit11_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function dT_value_Callback(hObject, eventdata, handles)
% hObject    handle to dT_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dT_value as text
%        str2double(get(hObject,'String')) returns contents of dT_value as a double
dT = str2double(get(hObject, 'String'));
if isnan(dT)
    set(hObject, 'String', 0);
    errordlg('Input must be a number','Error');
end

% Save the new Ryy_value value
handles.ukfdata.dT = dT;
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function dT_value_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dT_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function StateFcn_Callback(hObject, eventdata, handles)
% hObject    handle to StateFcn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of StateFcn as text
%        str2double(get(hObject,'String')) returns contents of StateFcn as a double
StateFcn = get(hObject, 'String');
if isempty(StateFcn)
    set(hObject, 'String', 0);
    errordlg('Input must be a number','Error');
end

% Save the new Ryy_value value
handles.ukfdata.StateFcn = StateFcn;
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function StateFcn_CreateFcn(hObject, eventdata, handles)
% hObject    handle to StateFcn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function MeasFcn_Callback(hObject, eventdata, handles)
% hObject    handle to MeasFcn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of MeasFcn as text
%        str2double(get(hObject,'String')) returns contents of MeasFcn as a double


% --- Executes during object creation, after setting all properties.
function MeasFcn_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MeasFcn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in Cfg_popup.
function Cfg_popup_Callback(hObject, eventdata, handles)
% hObject    handle to Cfg_popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns Cfg_popup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Cfg_popup


% --- Executes during object creation, after setting all properties.
function Cfg_popup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Cfg_popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function alphaSlider_Callback(hObject, eventdata, handles)
% hObject    handle to alphaSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
handles.ukfdata.alpha = get(hObject,'Value');
set(handles.alphaEdit,'string',num2str(handles.ukfdata.alpha));
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function alphaSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to alphaSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function bethaSlider_Callback(hObject, eventdata, handles)
% hObject    handle to bethaSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
handles.ukfdata.betha = get(hObject,'Value');
set(handles.bethaEdit,'string',num2str(handles.ukfdata.betha));
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function bethaSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bethaSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function kappaSlider_Callback(hObject, eventdata, handles)
% hObject    handle to kappaSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
handles.ukfdata.kappa = get(hObject,'Value');
set(handles.kappaEdit,'string',num2str(handles.ukfdata.kappa));
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function kappaSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to kappaSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on key press with focus on alphaSlider and none of its controls.
function alphaSlider_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to alphaSlider (see GCBO)
% eventdata  structure with the following fields (see UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)


function alphaEdit_Callback(hObject, eventdata, handles)
% hObject    handle to alphaEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of alphaEdit as text
%        str2double(get(hObject,'String')) returns contents of alphaEdit as a double
handles.ukfdata.alpha = str2double(get(hObject,'String'));
set(handles.alphaSlider,'Value',handles.ukfdata.alpha);
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function alphaEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to alphaEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function bethaEdit_Callback(hObject, eventdata, handles)
% hObject    handle to bethaEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of bethaEdit as text
%        str2double(get(hObject,'String')) returns contents of bethaEdit as a double
handles.ukfdata.betha = str2double(get(hObject,'String'));
set(handles.bethaSlider,'Value',handles.ukfdata.betha);
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function bethaEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bethaEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function kappaEdit_Callback(hObject, eventdata, handles)
% hObject    handle to kappaEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of kappaEdit as text
%        str2double(get(hObject,'String')) returns contents of kappaEdit as a double
handles.ukfdata.kappa = str2double(get(hObject,'String'));
set(handles.kappaSlider,'Value',handles.ukfdata.kappa);
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function kappaEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to kappaEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function stateMeaning_Callback(hObject, eventdata, handles)
% hObject    handle to stateMeaning (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of stateMeaning as text
%        str2double(get(hObject,'String')) returns contents of stateMeaning as a double


% --- Executes during object creation, after setting all properties.
function stateMeaning_CreateFcn(hObject, eventdata, handles)
% hObject    handle to stateMeaning (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupCfg.
function popupCfg_Callback(hObject, eventdata, handles)
% hObject    handle to popupCfg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupCfg contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupCfg
contents = cellstr(get(hObject,'String'));
handles.ukfdata.cfgid = str2double(contents{get(hObject,'Value')});
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function popupCfg_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupCfg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
