function varargout = untitled(varargin)

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @untitled_OpeningFcn, ...
                   'gui_OutputFcn',  @untitled_OutputFcn, ...
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


% --- Executes just before untitled is made visible.
function untitled_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);
set(handles.eq1, 'enable','off')
set(handles.e1, 'enable','off')
set(handles.eq2, 'enable','off')
set(handles.e2, 'enable','off')
set(handles.eq3, 'enable','off')
set(handles.e3, 'enable','off')
set(handles.eq4, 'enable','off')
set(handles.e4, 'enable','off')
set(handles.eq5, 'enable','off')
set(handles.e5, 'enable','off')
set(handles.GaussElimination, 'enable','off')
set(handles.GaussJordan, 'enable','off')
set(handles.GaussSeidel, 'enable','off')
set(handles.LU, 'enable','off')
set(handles.lhtext, 'enable','off')
set(handles.lhs, 'enable','off')
set(handles.filename, 'enable','off')
set(handles.name, 'enable','off')
set(handles.N, 'enable','off')
set(handles.enter, 'enable','off')
set(handles.enterfile, 'enable','off')
set(handles.dim, 'enable','off')
set(handles.tol, 'enable','off')
set(handles.dummy, 'enable','off')
set(handles.m, 'enable','off')
set(handles.t, 'enable','off')
set(handles.guesses, 'enable','off')
set(handles.dummy2, 'enable','off')
% UIWAIT makes untitled wait for user response (see UIRESUME)
% uiwait(handles.figure1);

function varargout = untitled_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;

function enter_Callback(hObject, eventdata, handles)
    global num 
    if num<=1
        msgbox('number of equations must be greater than 1');
        return;
    end
    set(handles.guesses, 'enable','on')
    set(handles.dummy2, 'enable','on')
    set(handles.GaussElimination, 'enable','on')
    set(handles.GaussJordan, 'enable','on')
    set(handles.GaussSeidel, 'enable','on')
    set(handles.LU, 'enable','on')
    set(handles.lhtext, 'enable','on')
    set(handles.lhs, 'enable','on')
    global tol
    global iterations;
    mynum = get(handles.N,'string');
    num = str2double(mynum);
    tol = get(handles.tol,'string');
    tol = str2double(tol);
    iterations = get(handles.m,'string');
    iterations = str2double(iterations);
    
    if num > 1
        set(handles.eq1, 'enable','on')
        set(handles.e1, 'enable','on')
    end
    if num >= 2
        set(handles.eq2, 'enable','on')
        set(handles.e2, 'enable','on')
    end
    if num >= 3
        set(handles.eq3, 'enable','on')
        set(handles.e3, 'enable','on')
        
    end
    if num >= 4
        set(handles.eq4, 'enable','on')
        set(handles.e4, 'enable','on')
        
    end
    if num >= 5
        set(handles.eq5, 'enable','on')
        set(handles.e5, 'enable','on')
    end
    

% --- Executes on button press in GaussElimination.
function GaussElimination_Callback(hObject, eventdata, handles)
    array = readvalues(hObject, eventdata, handles);
    if array(1)=='*'
        return;
    end
    output = Gauss_Elimination(array);
    if output(1)=='*'
        return;
    end
    write_in_gui(output,hObject,eventdata,handles);

    
    
function x = Gauss_Elimination(array)
    tic
    x=[];
    A= array;
    [m,n] = size(A);
    for j = 1:m-1
        for i= j+1:m
            A(i,:)=A(i,:)- A(j,:)*(A(i,j)/A(j,j));
        end  
    end 
    if A(m,n-1) == 0 && A(m,n)==0
        msgbox('The system has an infinite number of solutions');
        x(1)='*';
        return;
    end
    for i=1:m
        if A(i,i)==0
            msgbox('The system has no solutions');
            x(1)='*';
            return;
        end
    end

    x = zeros(1,m);
    for s =m:-1:1
        c=0;
        for k=2:m
            c=c+A(s,k)*x(k);
        end
        x(s) = (A(s,n)-c)/A(s,s);
    end
    time = num2str(toc);
    string1=['Time taken: ' time];
    msgbox(string1);
    x=x';
    writetable(array2table(x),'outputGauss.txt')
    

function x = Gauss_jordan(array)
    tic
    A= array ;
    [m,n] = size(A);
    x = zeros(1,m);
    for j = 1:m-1
        for i= j+1:m
            A(i,:)=A(i,:) - A (j,:)*(A(i,j)/A(j,j));
        end  
    end 
    for j=m:-1:2
        for i=j-1:-1:1
             A(i,:)=A(i,:) - A (j,:)*(A(i,j)/A(j,j));
        end
    end
    if A(m,n-1) == 0 && A(m,n)==0
        msgbox('The system has an infinite number of solutions');
        x(1)='*';
        return;
    end
    for i=1:m
        if A(i,i)==0
            msgbox('The system has no solutions');
            x(1)='*';
            return;
        end
    end

    for s=1:m
        A(s,:)=A(s,:)/A(s,s);
        x(s) = A(s,n);
    end
    time = num2str(toc);
    string1=['Time taken: ' time];
    msgbox(string1);
    x=x';
    writetable(array2table(x),'outputGaussJordan.txt');
            


function GaussJordan_Callback(hObject, eventdata, handles)
    global num;
    array = readvalues(hObject, eventdata, handles);
    if array(1)=='*'
        return;
    end
    output=Gauss_jordan(array);
    if output(1)=='*'
        return;
    end
    write_in_gui(output,hObject,eventdata,handles);

    function GaussSeidel_Callback(hObject, eventdata, handles)
    global num;
    array = readvalues(hObject, eventdata, handles);
    if array(1)=='*'
        return;
    end
    arr = get(handles.guesses,'string');
        t = regexp(regexprep(cellstr(arr), '^\s+', ''), '\s+', 'split');
        initial_guesses = str2double(vertcat(t{:}));
    L= length(initial_guesses);
    if L>num || L<num
        msgbox('number of initial guesses must agree with number of equations','error','Error');
        return;
    end
    output = Gauss_Seidel(array,initial_guesses);
    if output(1)=='*'
        return;
    end
    write_in_gui(output,hObject,eventdata,handles);
    
function x = checknum(array)
    global num
    L= length(array);
    x=0;
    if L>num || L<num
        x=1;
    end
   
    
    
function x1 = Gauss_Seidel(array, initial_g)
 tic
 global tol
 global iterations
 A = array;
 n = length(A)-1;
 a=A(:,1:n);
 R = sum(abs(a),2);
 D = abs(diag(a));
 W=R-D;
 check = D>=W;
 DD = all(check);
 if DD == 0
     msgbox('gauss seidel cannot be performed, matrix is not diagonally dominant', 'Error','error');
     x1(1)='*';
     return;
 end
 for i=1:n
     temp = array(i,i);
     if temp==0
         msgbox('gauss seidel cannot be performed diagonal contains a zero', 'Error','error');
         x1(1)='*';
         return;
     end
 end
 if det(a)==0
    msgbox('gauss seidel cannot be performed', 'Error','error');
    x1(1)='*';
    return;
 end
 
 x1 = initial_g;
 tolstring = isnan(tol);
 mycheck=all(tolstring);
 tempstring = isnan(iterations);
 mycheck2=all(tempstring);
 if mycheck==1
    tol = 0.00001;
 end
 m=50;
 if mycheck2==1
    m = 50;
 else
     m=iterations;
 end
 k = 1;
 X = zeros(n,m);
 errors = zeros(n,m);
 while  k <= m
    err = 0;
    for i = 1 : n 
       s = 0;
       for j = 1 : n 
          s = s-A(i,j)*x1(j);
       end
       s = (s+A(i,n+1))/A(i,i);
       if abs(s) > err 
         err  = abs(s);
       end
       errors(i,k)=err;
       x1(i) = x1(i) + s;
       X(i,k)=x1(i);
    end
    if err <= tol 
      break;
    else
      k = k+1;
    end
 end 
 time = num2str(toc);
 string1=['Time taken: ' time];
 string2=['  Number of iterations: ' num2str(k)];
 msgbox([string1 string2]);
 for i=1:n
     temparr = X(i,1:k);
     figure;
     n=linspace(0,k,k);
     plot(n,temparr);
     
 end
 output=x1;
 myerrors = errors(:,1:k);
 dummy = X(:,1:k);
 global num;
 N=num*2;
 table=zeros(k,N);
 i=1;
 v=1;
 while i<=N
     temp1 = myerrors(v,:)';
     temp2 = dummy(v,:)';
     table(:,i)=temp2;
     table(:,i+1)=temp1;
     i=i+2;
     v=v+1;
 end
 fig = figure;
 uit = uitable(fig,'Data',table);
 
 writetable(array2table(output'),'outputGaussSeidel.txt');
 
function LU_Callback(hObject, eventdata, handles)
    global num;
    array = readvalues(hObject, eventdata, handles);
    if array(1)=='*'
        return;
    end
    output = LUdec(array);
    if output(1)=='*'
        return;
    end
    write_in_gui(output,hObject,eventdata,handles);

function x= LUdec(array)
    tic
    A= array;
    [ccc,n] = size(A);
    A= array(:,1:n-1);
    B = array(:,n);
    [m,n]=size(A);
    if (m ~= n )
        msgbox( 'LU error: Matrix must be square' );
        x(1)='*';
        return; 
    end;
  L=zeros(m,m);
  U=zeros(m,m);
  for i=1:m
      % Finding L
      for k=1:i-1
          L(i,k)=A(i,k);
          for j=1:k-1
            L(i,k)= L(i,k)-L(i,j)*U(j,k);
          end
          L(i,k) = L(i,k)/U(k,k);
      end
      % Finding U
      for k=i:m
          U(i,k) = A(i,k);
          for j=1:i-1
          U(i,k)= U(i,k)-L(i,j)*U(j,k);
          end
      end
  end
  for i=1:m
      L(i,i)=1;   
  end
 
  y=zeros(m,1); % initiation for y 
  y(1)=B(1)/L(1,1);
  for i=2:m
      y(i)=-L(i,1)*y(1);
      for k=2:i-1
           y(i)=y(i)-L(i,k)*y(k);
      end;
          y(i)=(B(i)+y(i))/L(i,i);
   end;
    % Now we use this y to solve Ux = y
    x=zeros(m,1);
    if y(m)==0 && U(m,m)==0
        msgbox( 'LU error: system has infinte number of solutions' );
        x(1)='*';
        return;
    end
    for i=1:m
      temp = U(i,i);
      if temp==0
          msgbox( 'LU error: system has no solution' );
          x(1)='*';
          return;
      end
      
    end
    x(m)=y(m)/U(m,m);
    
    i=m-1;
    q=0;
    while  (i~= 0)
      x(i)=-U(i,m)*x(m);
       q=i+1;
          while (q~=m)
              x(i)=x(i)-U(i,q)*x(q);
              q=q+1;
          end;
        x(i)=(y(i)+x(i))/U(i,i);
        i=i-1;
    end
    time = num2str(toc);
    string1=['Time taken: ' time];
    msgbox(string1);
    output = x;
    writetable(array2table(output),'outputLU.txt');

function listbox2_Callback(hObject, eventdata, handles)
    %GUI deets%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    set(handles.tol, 'enable','on')
    set(handles.dummy, 'enable','on')
    set(handles.m, 'enable','on')
    set(handles.t, 'enable','on')
    set(handles.eq1, 'enable','off')
    set(handles.e1, 'enable','off')
    set(handles.eq2, 'enable','off')
    set(handles.e2, 'enable','off')
    set(handles.eq3, 'enable','off')
    set(handles.e3, 'enable','off')
    set(handles.eq4, 'enable','off')
    set(handles.e4, 'enable','off')
    set(handles.eq5, 'enable','off')
    set(handles.e5, 'enable','off')
    set(handles.GaussElimination, 'enable','off')
    set(handles.GaussJordan, 'enable','off')
    set(handles.GaussSeidel, 'enable','off')
    set(handles.LU, 'enable','off')
    set(handles.lhtext, 'enable','off')
    set(handles.lhs, 'enable','off')
    set(handles.filename, 'enable','off')
    set(handles.name, 'enable','off')
    set(handles.N, 'enable','off')
    set(handles.enter, 'enable','off')
    set(handles.enterfile, 'enable','off')
    set(handles.dim, 'enable','off')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    contents=get(handles.listbox2,'Value');
    
    if contents == 2
    set(handles.filename, 'enable','on')
    set(handles.name, 'enable','on')
    set(handles.enterfile, 'enable','on')
    end
    if contents == 1
    set(handles.N, 'enable','on')
    set(handles.enter, 'enable','on')
    set(handles.dim, 'enable','on')
    end
    



function enterfile_Callback(hObject, eventdata, handles)
    global num;
    global tol
    global iterations;
    tol = get(handles.tol,'string');
    tol = str2double(tol);
    iterations = get(handles.m,'string');
    iterations = str2double(iterations);
    filename = get(handles.name,'string');
    fileID = fopen(filename,'r');
    if fileID == -1
        warndlg('ERROR OPENING THE FILE');
    end
    if fileID ~= -1
    num = fgets(fileID);
    num = str2num(num);
    method = fgets(fileID);
    method = method(1:end-2);
    len =length(method);
    if len == 0
        msgbox('please specify the numerical method you want to use');
        return;
    end
    array =zeros(num,num+1);
    syms a b c d e
    
    i = 0;
    for z=1:num 
        arr = fgets(fileID);
        i = i + 1;
        t = regexp(regexprep(cellstr(arr), '^\s+', ''), '\s+', 'split');
        temparr = str2double(vertcat(t{:}));
        numnum = length(temparr);
        array(i,num+1) = temparr(numnum)*-1;
        if num >=1
            ca = coeffs(arr,a);
            temp = 0;
            if(length(ca)>1)
                temp = ca(2);
            end
            array(i,1) = temp;
        end

        if num >=2
            cb = coeffs(arr,b);
            temp = 0;
            if(length(cb)>1)
                temp = cb(2);
            end
            array(i,2) = temp;
        end
        if num >=3
            cc = coeffs(arr,c);
            temp = 0;
            if(length(cc)>1)
                temp = cc(2);
            end
            array(i,3) = temp;
        end
        if num >=4
            cd = coeffs(arr,d);
            array(i,4) = cd(2);
        end
        if num >=5
            ce = coeffs(arr,e);
            array(i,5) = ce(2);
        end 
    end
    g1=strcmpi(method, 'Gaussian-elimination');
    g2=strcmpi(method, 'Gaussian-jordan');
    g3=strcmpi(method, 'Gaussian-seidel');
    g4=strcmpi(method, 'LU-decomposition');
    [m,n] = size(array);
    if m~=n-1
        warndlg('MATRIX IS NOT SQUARE');
    end
    if m == n-1
    if g4 == 1
       output =  LUdec(array);
    end
    if  g3 == 1
        initialguesses = fgets(fileID);
        t= regexp(regexprep(cellstr(initialguesses), '^\s+', ''), '\s+', 'split');
        initialguesses = str2double(vertcat(t{:}));
        len = length(initialguesses);
        if len<num || len>num
            msgbox('Please make sure dimensions of initial guesses match number of equations');
            return
        end
        output = Gauss_Seidel(array, initialguesses);
    end
    if  g2 == 1
        output = Gauss_jordan(array);
    end
    if  g1 == 1
        output = Gauss_Elimination(array);
    end
        if output(1)=='*'
            return;
        end
        write_in_gui(output,hObject,eventdata,handles)
    end
    end
function ret = write_in_gui(output,hObject,eventdata,handles)
    global num;
        op=[];
    for i=1:num
        new = sprintf('\nVariable%d = %f\n\n',i, output(i));
        if i == 1
            op = new;
        end
        if i~=1
            op = [op new];
        end
    end
    set(handles.outputbox,'String',op);
   
    
function array = readvalues(hObject,eventdata,handles)
    global num;
    array =zeros(num,num+1);
    lhs_arr = get(handles.lhs,'string');
    t = regexp(regexprep(cellstr(lhs_arr), '^\s+', ''), '\s+', 'split');
    lhs_arr = str2double(vertcat(t{:}));
    x = checknum(lhs_arr);
    if x==1
       msgbox('LHS array dimensions do not match','error','Error');
       array(1)='*';
       return;
    end
    if num >= 1
        arr1 = get(handles.eq1,'string');
        t = regexp(regexprep(cellstr(arr1), '^\s+', ''), '\s+', 'split');
        arr1 = str2double(vertcat(t{:}));
        x = checknum(arr1);
        if x==1
            msgbox('First equation dimensions do not match','error','Error');
            array(1)='*';
            return;
        end
        array = arr1;
    end 
    if num >= 2
        arr2 = get(handles.eq2,'string');
        t = regexp(regexprep(cellstr(arr2), '^\s+', ''), '\s+', 'split');
        arr2 = str2double(vertcat(t{:}));
        x = checknum(arr2);
        if x==1
            msgbox('Second equation dimensions do not match','error','Error');
            array(1)='*';
            return;
        end
        array = [array; arr2];
    end 
    if num >= 3
        arr3 = get(handles.eq3,'string');
        t = regexp(regexprep(cellstr(arr3), '^\s+', ''), '\s+', 'split');
        arr3 = str2double(vertcat(t{:}));
        x = checknum(arr3);
        if x==1
            msgbox('Third equation dimensions do not match','error','Error')
            array(1)='*';
            return;
        end
        array = [array; arr3];
    end 
    if num >= 4
        arr4 = get(handles.eq4,'string');
        t = regexp(regexprep(cellstr(arr4), '^\s+', ''), '\s+', 'split');
        arr4 = str2double(vertcat(t{:}));
        x = checknum(arr4);
        if x==1
            msgbox('Fourth equation dimensions do not match','error','Error');
            array(1)='*';
            return;
        end
        array = [array; arr4];
    end 
    if num >= 5
        arr5 = get(handles.eq5,'string');
        t = regexp(regexprep(cellstr(arr5), '^\s+', ''), '\s+', 'split');
        arr5 = str2double(vertcat(t{:}));
        x = checknum(arr5);
        if x==1
            msgbox('Fifth equation dimensions do not match','error','Error');
            array(1)='*';
            return;
        end
        array = [array;arr5];
    end
    b = lhs_arr';
    array = [array b];

function tol_Callback(hObject, eventdata, handles)
function tol_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function m_Callback(hObject, eventdata, handles)
function m_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function guesses_Callback(hObject, eventdata, handles)
function guesses_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function listbox2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function N_Callback(hObject, eventdata, handles)
function N_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function eq1_Callback(hObject, eventdata, handles)
function eq1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function eq2_Callback(hObject, eventdata, handles)
function eq2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function eq3_Callback(hObject, eventdata, handles)
function eq3_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function eq4_Callback(hObject, eventdata, handles)
function eq4_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function eq5_Callback(hObject, eventdata, handles)
function eq5_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function name_Callback(hObject, eventdata, handles)
function name_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function lhs_Callback(hObject, eventdata, handles)

function lhs_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
