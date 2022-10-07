function varargout = MTEX2GmshGUI(varargin)
% MTEX2GMSHGUI MATLAB code for MTEX2GmshGUI.fig
%      MTEX2GMSHGUI, by itself, creates a new MTEX2GMSHGUI or raises the existing
%      singleton*.
%
%      H = MTEX2GMSHGUI returns the handle to a new MTEX2GMSHGUI or the handle to
%      the existing singleton*.
%
%      MTEX2GMSHGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MTEX2GMSHGUI.M with the given input arguments.
%
%      MTEX2GMSHGUI('Property','Value',...) creates a new MTEX2GMSHGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before MTEX2GmshGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to MTEX2GmshGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help MTEX2GmshGUI

% Last Modified by GUIDE v2.5 07-Oct-2022 11:19:01

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @MTEX2GmshGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @MTEX2GmshGUI_OutputFcn, ...
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



% --- Executes just before MTEX2GmshGUI is made visible.
function MTEX2GmshGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to MTEX2GmshGUI (see VARARGIN)

% Choose default command line output for MTEX2GmshGUI
handles.output = hObject;
if nargin>3 && isa(varargin{1}, 'gmshGeo')
    handles.G=varargin{1};
else
    handles.G=[];
end

% Update handles structure
guidata(hObject, handles);

% Import gmshGeo variables from workspace to GUI
evalin('base','s_tmp=whos;');
my_variables = evalin('base','{s_tmp(strcmp({s_tmp.class}, ''gmshGeo'')).name}');
if isempty(my_variables)
    set(handles.varname,'String',{'-No gmshGeo object found in the workspace-'});
    set(handles.plot,'enable','off');
    set(handles.exportgrainprops,'enable','off');
else
    set(handles.varname,'String',my_variables);
end
evalin('base','clear s_tmp');




% UIWAIT makes MTEX2GmshGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = MTEX2GmshGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function filepath_Callback(hObject, eventdata, handles)
% hObject    handle to filepath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of filepath as text
%        str2double(get(hObject,'String')) returns contents of filepath as a double


% --- Executes during object creation, after setting all properties.
function filepath_CreateFcn(hObject, eventdata, handles)
% hObject    handle to filepath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in browse.
function browse_Callback(hObject, eventdata, handles)
% hObject    handle to browse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
filter={'*.msh','Gmsh MSH  (*.msh)'; ...
        '*.*',  'Other mesh formats (*.*)'};
[filename, pathname]=uiputfile(filter,'Save mesh as');
set(handles.filepath,'string',[pathname filename])


% --- Executes on button press in medium.
function medium_Callback(hObject, eventdata, handles)
% hObject    handle to medium (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of medium
textfields={'mediumx' 'mediumy' 'mediumz' 'MediumElementSize'};
if strcmp(get(handles.mediumx,'enable'),'off')
    set_to='on';
else
    set_to='off';
end
for i=1:4
    set(handles.(textfields{i}),'enable',set_to);
end




function mediumx_Callback(hObject, eventdata, handles)
% hObject    handle to mediumx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of mediumx as text
%        str2double(get(hObject,'String')) returns contents of mediumx as a double


% --- Executes during object creation, after setting all properties.
function mediumx_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mediumx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function mediumy_Callback(hObject, eventdata, handles)
% hObject    handle to mediumy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of mediumy as text
%        str2double(get(hObject,'String')) returns contents of mediumy as a double


% --- Executes during object creation, after setting all properties.
function mediumy_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mediumy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function mediumz_Callback(hObject, eventdata, handles)
% hObject    handle to mediumz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of mediumz as text
%        str2double(get(hObject,'String')) returns contents of mediumz as a double


% --- Executes during object creation, after setting all properties.
function mediumz_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mediumz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function MediumElementSize_Callback(hObject, eventdata, handles)
% hObject    handle to MediumElementSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of MediumElementSize as text
%        str2double(get(hObject,'String')) returns contents of MediumElementSize as a double


% --- Executes during object creation, after setting all properties.
function MediumElementSize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MediumElementSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function cmdbak_Callback(hObject, eventdata, handles)
% hObject    handle to cmdbak (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of cmdbak as text
%        str2double(get(hObject,'String')) returns contents of cmdbak as a double


% --- Executes during object creation, after setting all properties.
function cmdbak_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cmdbak (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in run.
function run_Callback(hObject, eventdata, handles)
% hObject    handle to run (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
cmd=updateCmdLine(handles);
G=loadG(hObject, eventdata, handles);
eval(cmd)



function ElementOrder_Callback(hObject, eventdata, handles)
% hObject    handle to ElementOrder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ElementOrder as text
%        str2double(get(hObject,'String')) returns contents of ElementOrder as a double


% --- Executes during object creation, after setting all properties.
function ElementOrder_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ElementOrder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ElementSize_Callback(hObject, eventdata, handles)
% hObject    handle to ElementSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ElementSize as text
%        str2double(get(hObject,'String')) returns contents of ElementSize as a double


% --- Executes during object creation, after setting all properties.
function ElementSize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ElementSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in usegradient.
function usegradient_Callback(hObject, eventdata, handles)
% hObject    handle to usegradient (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of usegradient
if strcmp(get(handles.gradient,'enable'),'off')
    set(handles.gradient,'enable','on');
    set(handles.usecurvature,'Value',0);
    set(handles.Curvature,'string','0');
    set(handles.Curvature,'enable','off');
else
    set(handles.gradient,'enable','off');
    set(handles.gradient,'string','0');
end
    


function slope_Callback(hObject, eventdata, handles)
% hObject    handle to usegradient (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of usegradient as text
%        str2double(get(hObject,'String')) returns contents of usegradient as a double


% --- Executes during object creation, after setting all properties.
function usegradient_CreateFcn(hObject, eventdata, handles)
% hObject    handle to usegradient (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function gradient_CreateFcn(hObject, eventdata, handles)
% hObject    handle to usegradient (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in usecurvature.
function usecurvature_Callback(hObject, eventdata, handles)
% hObject    handle to usecurvature (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of usecurvature
if strcmp(get(handles.Curvature,'enable'),'off')
    set(handles.Curvature,'enable','on');
    set(handles.Curvature,'string','4');
    set(handles.usegradient,'Value',0);
    set(handles.usegradient,'enable','off');
    set(handles.gradient,'string','0');
    set(handles.gradient,'enable','off');    
else
    set(handles.Curvature,'enable','off');
    set(handles.Curvature,'string','0');
    set(handles.usegradient,'enable','on');
end


% --- Executes during object creation, after setting all properties.
function Curvature_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Curvature (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in exportgrainprops.
function exportgrainprops_Callback(hObject, eventdata, handles)
% hObject    handle to exportgrainprops (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
filter={'*.csv','Comma-separated values  (*.csv)'; ...
        '*.txt',  'ASCII Text file (*.txt)'};
[filename, pathname]=uiputfile(filter,'Save mesh as');
cols=cell(3,1);
colnames={'Phase', 'PhaseID', 'Euler', 'Rodrigues', 'GrainID'};
k=0;
for i=1:3
    colname=sprintf('Column%i',i);
    id_selected=get(handles.(colname),'value');
    if id_selected>5
        break
    end
    k=k+1;
    colval =colnames{id_selected};
    cols{i} = sprintf('''%s''',strrep(colval,' ',''));
end
list=join(cols(1:k),',');
G=loadG(hObject, eventdata, handles);
cmd=sprintf('G.exportGrainProps(%s, %s)',[pathname filename], list{1});
    



function Curvature_Callback(hObject, eventdata, handles)
% hObject    handle to Curvature (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Curvature as text
%        str2double(get(hObject,'String')) returns contents of Curvature as a double

function gradient_Callback(hObject, eventdata, handles)
% hObject    handle to Curvature (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Curvature as text
%        str2double(get(hObject,'String')) returns contents of Curvature as a double



% --- Executes during object creation, after setting all properties.
function Curv_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Curvature (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




function verbosity_Callback(hObject, eventdata, handles)
% hObject    handle to verbosity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of verbosity as text
%        str2double(get(hObject,'String')) returns contents of verbosity as a double


% --- Executes during object creation, after setting all properties.
function verbosity_CreateFcn(hObject, eventdata, handles)
% hObject    handle to verbosity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function partition_Callback(hObject, eventdata, handles)
% hObject    handle to partition (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of partition as text
%        str2double(get(hObject,'String')) returns contents of partition as a double


% --- Executes during object creation, after setting all properties.
function partition_CreateFcn(hObject, eventdata, handles)
% hObject    handle to partition (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Layers_Callback(hObject, eventdata, handles)
% hObject    handle to Layers (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Layers as text
%        str2double(get(hObject,'String')) returns contents of Layers as a double


% --- Executes during object creation, after setting all properties.
function Layers_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Layers (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function cmd=updateCmdLine(handles)
    filepath=get(handles.filepath,'String');
    options_fields={'ElementSize', 'Thickness', 'gradient', 'ElementType', 'ElementOrder', 'Periodic', 'Curvature', 'grainPrefix', 'Medium', 'verbosity', 'partition', 'MediumElementSize', 'Layers', 'LocalSize'};
    option_default={[], [], 0, 'Wedge',   1   , 'None', 0, 'Grain_', [0, 0, 0], 4, 1, [], 1, []};
    nopt=length(options_fields);
    cmds=cell(nopt,1);
    k=0;
    for i=1:nopt
        field=options_fields{i};
        if strcmp(field,'ElementType')
            value=get(get(handles.(field),'SelectedObject'),'Tag');
        elseif strcmp(field,'Periodic')
            if get(handles.Xperiodic,'Value')
                if get(handles.Yperiodic,'Value')
                    value='Both';
                else
                    value='X';
                end
            elseif get(handles.Yperiodic,'Value')
                value='Y';
            else
                value='None';
            end
        elseif strcmp(field,'Medium')
            if get(handles.medium,'Value')
                value=[str2double(get(handles.mediumx,'String')) str2double(get(handles.mediumy,'String')) str2double(get(handles.mediumz,'String'))];
            else
                value=option_default{i};
            end
        elseif strcmp(field,'MediumElementSize')
            if ~isempty(get(handles.MediumElementSize,'string')) && get(handles.medium,'Value')
                value=str2double(get(handles.MediumElementSize,'String'));
            else
                value=option_default{i};
            end
        elseif strcmp(field,'LocalSize')
            if get(handles.LocalSize,'Value')
                value=get(handles.LocalSizesValues, 'data');
            else
                value=[];
            end
        else
            value=get(handles.(field),'String');
        end
        if isnumeric(option_default{i}) && ischar(value)
            value=str2double(value);    % If the default value is numeric, convert the new value
        end
        
        % Append the Name-value pair arguments if different from default
        if isa(value,'cell') || ( ~isequal(value, option_default{i}) && ~any(isnan(value)) )
            k=k+1;
            if isa(value,'cell')    % Use LocalSize function
                m=cell2mat(value);
                s1=sprintf('%i ',m(:,1));
                s2=sprintf('%g ',m(:,2));
                s3=sprintf('%g ',m(:,3));
                value=sprintf('LocalSize([%s],[%s],[%s])',s1,s2,s3);
            elseif isnumeric(value)
                if length(value)==1
                    value=sprintf('%d',value);
                else
                    value=sprintf('[%d %d %d]',value);
                end 
            else
                value=sprintf('''%s''',value);
            end
            cmds{k}=[sprintf('''%s''',options_fields{i}) ',' value];
        end
    end
    if k % Do this only if arguments have to be appended
        cmds=join(string(cmds(1:k)),',');
        cmd= [sprintf('G.mesh(''%s''', filepath) ',' cmds{1} ')'];
    else
        cmd= [sprintf('G.mesh(''%s'')', filepath)];
    end

function grainPrefix_Callback(hObject, eventdata, handles)
% hObject    handle to grainPrefix (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of grainPrefix as text
%        str2double(get(hObject,'String')) returns contents of grainPrefix as a double


% --- Executes during object creation, after setting all properties.
function grainPrefix_CreateFcn(hObject, eventdata, handles)
% hObject    handle to grainPrefix (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when selected object is changed in ElementType.
function ElementType_SelectionChangedFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in ElementType 
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function Thickness_Callback(hObject, eventdata, handles)
% hObject    handle to Thickness (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Thickness as text
%        str2double(get(hObject,'String')) returns contents of Thickness as a double


% --- Executes during object creation, after setting all properties.
function Thickness_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Thickness (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in Periodic.
function Periodic_Callback(hObject, eventdata, handles)
% hObject    handle to Periodic (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns Periodic contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Periodic


function Xperiodic_Callback(hObject, eventdata, handles)
% hObject    handle to usegradient (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of usegradient

function Yperiodic_Callback(hObject, eventdata, handles)
% hObject    handle to usegradient (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of usegradient


% --- Executes on button press in showcommand.
function showcommand_Callback(hObject, eventdata, handles)
% hObject    handle to showcommand (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
txt=updateCmdLine(handles);
answer = questdlg(txt, 'Command line', 'Copy to clipboard','Close','Close');
% Handle response
if strcmp(answer,'Copy to clipboard')
    clipboard('copy',txt)
end


% --- Executes on selection change in Column1.
function Column1_Callback(hObject, eventdata, handles)
% hObject    handle to Column1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns Column1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Column1


% --- Executes during object creation, after setting all properties.
function Column1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Column1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in Column2.
function Column2_Callback(hObject, eventdata, handles)
% hObject    handle to Column2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns Column2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Column2


% --- Executes during object creation, after setting all properties.
function Column2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Column2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in Column3.
function Column3_Callback(hObject, eventdata, handles)
% hObject    handle to Column3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns Column3 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Column3


% --- Executes during object creation, after setting all properties.
function Column3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Column3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in LocalSizesValues.
function LocalSize_Callback(hObject, eventdata, handles)
% hObject    handle to LocalSizesValues (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of LocalSizesValues


% --- Executes on selection change in varname.
function varname_Callback(hObject, eventdata, handles)
% hObject    handle to varname (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns varname contents as cell array
%        contents{get(hObject,'Value')} returns selected item from varname


% --- Executes during object creation, after setting all properties.
function varname_CreateFcn(hObject, eventdata, handles)
% hObject    handle to varname (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in plot.
function plot_Callback(hObject, eventdata, handles)
% hObject    handle to plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA
G=loadG(hObject, eventdata, handles);
figure
plot(G);

function G=loadG(hObject, eventdata, handles)
varnames=get(handles.varname,'string');
G = evalin('base',varnames{get(handles.varname,'value')});
