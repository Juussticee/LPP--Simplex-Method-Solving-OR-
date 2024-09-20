    clc
clear all

% First Window: Read cost function, max/min, number of constraints
options.Resize = 'on';
title = 'Enter your LP data';
prompt = {'Cost function (c vector)',...
    'Max (enter 1) Min (enter 2)',...
    'Number of constraints'};
DefaultText = {'[ ]','',''};
a = inputdlg(prompt,title,1,DefaultText,options);

c = eval(a{1});
type = eval(a{2});
nConstraints = eval(a{3});
struct1 = struct('vari',{},'Type',{});

% Second window: Sense of each inequality constraint
title = '(In)equality Constraints Data';
for i = 1:nConstraints
    prompt = {['Specify (<=,>=,=) for constraint no. ', num2str(i)]};
    DefaultText = {''};
    a = inputdlg(prompt,title,1,DefaultText,options);
    struct1(1,i).Type = a{1};
end

% Third window: matrix A
title = 'Inequality constraint coefficients';
prompt = {'GIVE ME MATLAB MATRIX A the variables values (x1 x2 ...etc)'};
DefaultText = {'[]'}; % MATLAB MATRIX
a = inputdlg(prompt,title,1,DefaultText,options);
A = eval(a{1});

% Fourth window: vector b
title = 'Inequality (upper/lower) Bounds';
prompt = {'GIVE ME VECTOR b'};
DefaultText = {'[ ]'}; % VECTOR MATRIX
a = inputdlg(prompt,title,1,DefaultText,options);
b = eval(a{1});

% Slack, Surplus, Non-basic (activity)
slack = []; % Slack
Aslack = []; % Matrix for slacks
artificial = []; % Artificial or surplus
AArtificial = []; % Matrix for surplus
v_b = []; % Vector of bounds

for i = 1:nConstraints
    n = struct1(1,i).Type;
    
    switch n
        case '<='
            slack = [slack b(i)];
            Aslack(i,length(slack)) = 1;
            
        case '>='
            slack = [slack 0];
            Aslack(i,length(slack)) = -1;
            artificial =[artificial b(i)];
            AArtificial(i,length(artificial)) = 1;
            
        case '='
            artificial =[artificial b(i)];
            AArtificial(i,length(artificial)) = 1;
    end
    v_b = [v_b b(i)];
end
disp(size(A));
disp(size(Aslack));
disp(size(AArtificial));
% Adjust dimensions of Aslack
if size(A, 1) > size(Aslack, 1)
    Aslack = [Aslack; zeros(size(A, 1) - size(Aslack, 1), size(Aslack, 2))];
end

% Adjust dimensions of AArtificial
if size(A, 1) > size(AArtificial, 1)
    AArtificial = [AArtificial; zeros(size(A, 1) - size(AArtificial, 1), size(AArtificial, 2))];
end

% Concatenate matrices A, Aslack, and AArtificial
A = [A, Aslack, AArtificial];
vari = [];
vari_artificial = []; % Artificials (surplus)
vari_slack = []; % Slacks
vari_non = []; % Non-basic (activity)

for i = 1:size(A,2)
    struct1(1,i).vari = ['x',num2str(i)];
    vari = [vari, struct1(1,i).vari,' '];
    if i <= length(c)
        vari_non = [vari_non, struct1(1,i).vari,' ']; % Activity variables x1 , x2
    elseif i <= length(c) + length(slack)
        vari_slack = [vari_slack, struct1(1,i).vari,' ']; % Slack
    else
        vari_artificial = [vari_artificial, struct1(1,i).vari,' ']; % Artificial
    end
end

v_a = zeros(size(c));

x = [v_a,slack,artificial];
%if isempty(artificial)
 %   v_ar = [];
%else
 %   if type == 1 % Max
  %      v_ar = -1 * ones(1,length(artificial));
   % else
    %    v_ar = +1 * ones(1,length(artificial));
    %end
%end
if type == 2 %Min
    c=-c;
end 
c1=c;
for i=1:length(slack)
    c1=[c1,0];
end 
for i=1:length(artificial)
c1=[c1,1];
end 
Z1=[];
for i=1:length(slack);
    Z1=[Z1,0];
end 
for i=1:length(artificial);
    Z1=[Z1,1];
end 
% Initialize Cj,Zj

    Cj = ['Z  |', vari,'|  Solution  |' , 'Ratio' ];
    Ci = ['Z'];
    for i=1:length(artificial)
        Ci=[Ci ; "a"+num2str(i)];
    end 
    Cs=[];
    tabl = [];
    Vb = [];
    Q = v_b;
    A=[c1;A];
    for i = 1:length(Cj)
        tabl = [tabl; ' | '];
        struct2(1,i).value = Q(i);
        ind = find(x == Q(i));
        struct2(1,i).var_base = struct1(1,ind).vari;
        Vb = [Vb,struct2(1,i).var_base,' '];
        Cs = [Cs, Cj(ind)];
    end

Z = sum(Cs.*Q);
for i = 1:length(Cj)
    Zj(i) = sum(Cs'.*A(:,i));
end

vars_at_moment=[];
for i = 1:nConstraints
    if(length(struct2(1,i).var_base) == 2)
        vars_at_moment=[vars_at_moment;struct2(1,i).var_base,' '];
    else
        vars_at_moment=[vars_at_moment;struct2(1,i).var_base];
    end
end

fprintf('\n')

disp('==== LP in Standard Form =====')
disp(['Variables: ', vari])
disp(['     -Activity (non-basic) Variables: ', vari_non])
disp(['     -Slack Variables: ', vari_slack])
disp(['     -Artificial (surplus) Variables: ', vari_artificial])

disp('==== Iteration 0 =====')
disp([' Initializing my variables: ', vari])
disp(['     -Activity (non-basic) Variables: ', num2str(v_a)])
disp(['     -Slack Variables: ', num2str(slack)])
disp(['     -Artificial (surplus) Variables: ', num2str(v_ar)])
disp('===================')
disp(['Cj: ',num2str(Cj)])
disp([tabl,num2str(Cs'),tabl,vars_at_moment,tabl,num2str(Q'),tabl,num2str(A),tabl])
disp('===================')
disp(['Zj: ',num2str(Zj)])
disp(['Cj-Zj: ',num2str(Cj-Zj)])
disp(['Z: ',num2str(Z)])

iterNum = 1;
maxIterations = 10; % Define maximum number of iterations

while iterNum <= maxIterations
    % Entering variable
    if type == 1 % Max
        num = max(Cj - Zj);
        num = num(1);
        num1 = find(Cj - Zj == num);
        num1 = num1(1);
        V_enter = struct1(1,num1);
    else % Min
        num = min(Cj - Zj);
        num = num(1);
        num1 = find(Cj - Zj == num);
        num1 = num1(1);
        V_enter = struct1(1,num1);
    end
    
    b = A(:,num1);
    k = 0;
    d = 1e4;
    for i = 1:length(Q)
       if b(i) > 0
          div = Q(i)/b(i);
          if d > div
             d = div;
             k = i;
          end
       end
    end
    if k == 0
        disp('Solution is infinity');
        break;
    else
        num2 = k;
    end
    
    V_leave = struct2(1,num2).var_base;
    struct2(1,num2).var_base = struct1(1,num1).vari;
    pivot = A(num2,num1);
    Cs(num2) = Cj(num1);
    A(num2,:) = A(num2,:)./pivot;
    Q(num2) = Q(num2)./pivot;
    h = size(A,1);
    for i = 1:h
        if i ~= num2
           Q(i) = Q(i) - A(i,num1)*Q(num2);
           A(i,:) = A(i,:) - A(i,num1)*A(num2,:);
        end
    end
    
    Z = sum(Cs.*Q);
    for i = 1:size(A,2)
       Zj(i) = sum(Cs'.*A(:,i)); 
    end
    
    vars_at_moment=[];
    for i = 1:nConstraints
        if(length(struct2(1,i).var_base) == 2)
            vars_at_moment=[vars_at_moment;struct2(1,i).var_base,' '];
        else
            vars_at_moment=[vars_at_moment;struct2(1,i).var_base];
        end
    end
    
    disp(['==== Iteration', num2str(iterNum) ,'=====']);
    disp(['ENTER: ', V_enter.vari]);
    disp(['LEAVE: ', V_leave]);
    disp(['PIVOT: ', num2str(pivot)]);
    
    disp('===================')
    disp(['Cj: ',num2str(Cj)])
    disp([tabl,num2str(Cs'),tabl,vars_at_moment,tabl,num2str(Q'),tabl,num2str(A),tabl])
    disp('===================')
    disp(['Zj: ',num2str(Zj)])
    disp(['Cj-Zj: ',num2str(Cj-Zj)])
    disp(['Z: ',num2str(Z)])
    
    % Stopping criterion
    if type == 1 % Max
       temp = max(Cj - Zj);
       temp = temp(1);
       if temp <= 0
          break; 
       end
    else % Min
        temp = min(Cj - Zj);
        temp = temp(1);
       if temp >= 0
          break; 
       end 
    end
    
    iterNum = iterNum + 1;
end
