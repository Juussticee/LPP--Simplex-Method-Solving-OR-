clc;
clear all;
format short;
% First Window: Read cost function, max/min, number of constraints
options.Resize = 'on';%%so we can resize the GUI box that appears
title = 'Enter your LP data';
prompt = {'Cost function (c vector)',...
    'Max (enter 1) Min (enter 2)',...
    'Number of constraints'};
DefaultText = {'[ ]','',''};
a = inputdlg(prompt,title,1,DefaultText,options);%%inputdllg will open the dialogue box

c = eval(a{1});%%eval its like num2str or num2double but wider which mean it can accept any value whether string or array
type = eval(a{2});%%or even commands and that what makes it risky because some entered commands can risk the device,
nConstraints = eval(a{3});%%also it take more time to execute and compile then num2str and those. it's rec. for advanced codes
struct1 = struct('vari',{},'Type',{});

% Second window: Sense of each inequality constraint
title = '(In)equality Constraints Data';
for i = 1:nConstraints
    prompt = {['Specify (<=,>=,=) for constraint no. ', num2str(i)]};
    DefaultText = {''};%%how the text neterd is expected to be
    a = inputdlg(prompt,title,1,DefaultText,options);
    struct1(1,i).Type = a{1};%% structured array , it's like entities or classes where u can store many information
end %%and data in one , like if you want to assign a preson with a name and height etc... it'll beeasier to assign them
%all in a structured array specified for that person


% Third window: matrix A
title = 'Inequality constraint coefficients';
prompt = {'GIVE ME MATLAB MATRIX A the variables values (x1 x2 ...etc)'};
DefaultText = {'[]'}; % MATLAB MATRIX
a = inputdlg(prompt,title,1,DefaultText,options);
A = eval(a{1});
% Fourth window: vector b (answers)
title = 'Inequality (upper/lower) Bounds';
prompt = {'GIVE ME VECTOR b'};
DefaultText = {'[]'}; % VECTOR MATRIX
a = inputdlg(prompt,title,1,DefaultText,options);
b = eval(a{1});

% Slack, Surplus, Non-basic (activity)
slack = []; % Slack
Aslack = []; % Matrix for slacks
artificial = []; % Artificial or surplus
AArtificial = []; % Matrix for surplus
v_b = []; % Vector of bounds
vars=struct('Nb',{});
vars(1,1).Nb=['Z'];
for i = 1:nConstraints
    n = struct1(1,i).Type;
    switch n
        case '<='
            slack = [slack b(i)];%%Compute the values of the variables
            Aslack(i,length(slack)) = 1;
            vars(1,i+1)   .Nb=['s',num2str(i)];
            operation_type=1;%%Simplex
        case '>='
            slack = [slack b(i)];
            Aslack(i,length(slack)) =-1;
            artificial =[artificial b(i)];
            AArtificial(i,length(artificial)) = 1;
            vars(1,i+1).Nb=['a',num2str(i)];
            operation_type=2;%%Two-Phase-Method
            
        case '='
            artificial =[artificial b(i)];
            AArtificial(i,length(artificial)) = 1;
            vars(1,i+1).Nb=['a',num2str(i)];
            operation_type=2;%%Two-Phase-Method
    end
    v_b = [v_b b(i)];
end
disp([size(A),'size of x']);
disp([size(Aslack),'size of slacks']);
disp([size(AArtificial),'size of artificial']);
% Adjust dimensions of Aslack
if size(A, 1) > size(Aslack, 1)
    Aslack = [Aslack; zeros(size(A, 1) - size(Aslack, 1), size(Aslack, 2))];
end
% Adjust dimensions of AArtificial
if size(A, 1) > size(AArtificial, 1)
    AArtificial = [AArtificial; zeros(size(A, 1) - size(AArtificial, 1), size(AArtificial, 2))];
end%% so we can add them to A

% Concatenate matrices A, Aslack, and AArtificial
Aa=[];
A = [A, Aslack, AArtificial];
Aa=[A,transpose(v_b)];%%here we obtain the table we usually start working on
disp(Aa);
vari = [];
vari_artificial = []; % Artificials
vari_slack = []; % Slacks
vari_non = []; % Non-basic (activity)
non=[];
varis_slack=[];
varis_artificial=[];
for i = 1:size(A,2)%%we name the variab;les , x1 , x2 
    struct1(1,i).vari = ['x',num2str(i)];
    vari = [vari, struct1(1,i).vari,' '];
    if i <= length(c)
        vari_non = [vari_non, struct1(1,i).vari,' ']; % Activity variables x1 , x2
        non=[non,i];%%we assign each var is where (so we can find later)
    elseif i <= length(c) + length(slack)
        vari_slack = [vari_slack, struct1(1,i).vari,' ']; % Slack
        varis_slack=[varis_slack,i];
    else
        vari_artificial = [vari_artificial, struct1(1,i).vari,' ']; % Artificial
        varis_artificial=[varis_artificial,i];
    end
end
c1=c;
if type == 2 %Min
    c=-c;
end
ci=c;
for i=1:length(slack)
    c1=[c1,0];
    ci=[ci,0];
end 
ci=[ci,0]
for i=1:length(artificial)
    c1=[c1,1];
end
    varis=strsplit(vari);
    varis=varis(~cellfun('isempty', varis));%%doing this to switch them to cells itll be easier to get the length value like this 
 %%
 c2=[c1,0];%% just for backup , keeping c1 aside
 Aa=[c2;Aa];
 varis_non=strsplit(vari_non);
 varis_non=varis_non(~cellfun('isempty',varis_non));
 Aa2=[Aa];   %%NOte 2-phase method isnt complete yet 
 if operation_type==2%%if its 2 phase method 
 for i=1:length(varis_non)
     Aa2(1,i)=0;
 end
 for i=1:length(artificial)
     r_nb=find(1==Aa2(:,varis_artificial(1,i)));
     Aa2(1,:)=Aa2(1,:)-Aa2(r_nb(2,1),:);
 end
 end 
 for i=1:length(slack)
     disp(['S', num2str(i), '=', num2str(slack(1, i))]);
 end 
 for i=1:length(artificial)
     disp(['a', num2str(i), '=', num2str(slack(1, i))]);
 end 
 %%
 pivot=0;
 vars_char = char(vars.Nb);
 %vars_cell=cellstr(vars_char);
 ar_column=[];
 Aa3=[Aa2];
     if length(non)+length(slack)+1<length(Aa3);
         art_fnd=1;
     else art_fnd=0;
     end 
     negative_value1=min(Aa3(1,1:length(Aa3)-1))<0;
while art_fnd>0
       while negative_value1==true
         expt=length(artificial)+1;
         enter_var=find(min(Aa3(1,1:end-expt))==Aa3(1,:),1);
         leave_var=0;
         ratio1=[];
         for k = 2:height(Aa)
            if Aa3(k, end) / Aa3(k, enter_var) <= 0
                 ratio1 = [ratio1,0];
            else
                 ratio1 = [ratio1, Aa3(k, end) / Aa3(k, enter_var)];
            end
         end
          leave_var = find(ratio1 == min(ratio1),1)+1;
          pivot=Aa3(leave_var,enter_var);%%%%%
          for g=1:length(Aa3)
          Aa3(leave_var,g)=[Aa3(leave_var,g)/pivot];
          end 
          vars_char(leave_var,:)=char(varis(1,enter_var));
          pivot_row=leave_var;
          for t=1:height(Aa3)
              if t~=pivot_row
                 remainder=Aa3(t,enter_var);
                Aa3(t,:)=Aa3(t,:)-(remainder*Aa3(leave_var,:));
              end 
          end 
          delete=height(Aa3)-pivot_row
          del_col_num=length(Aa3)-1-pivot_row;
          Aa3(:, del_col_num) = []; 
          negative_value1=min(Aa3(1,1:length(Aa3)-1))<0;
        if length(non)+length(slack)+1<length(Aa3);
           art_fnd=1;
        else art_fnd=0;
        end
     end 
end  
if operation_type==2
for i=1:length(non)
    if Aa2(1,i)>0
        r_nb=find(1==Aa2(:,non(1,i)));
            Aa2(1,:)=Aa2(1,:)-Aa2(r_nb(1,1),:);
    end 
end
end 
%
%
%
     %Display the results
     if operation_type==2
     disp('---Optimal solution to 1st phase =');
     disp(['Z, max = ',Aa3(1,end)]);
     for i=2:height(Aa3)
     disp(['x',num2str(i-1)]); %  we do -1 because the height starts from 2
     disp(['=',num2str(Aa3(i,end))]);
     end
     end 
     %%@nd phase or simplex method 
     As1=[Aa3];
     for i=1:length(varis_non);
         As1(1,i)=ci(1,i);
     end
     %% begin simplex
     As2=[As1];
     negative_value=min(As2(1,1:end-1))<0;
     while negative_value==true
         negative_value=0;%%reset all their values to 0 cause sometimes they were doing errors and not reseting 
         min_value=0;
         enter_var=0;
         min_ratio_value=0;
         min_value=min(As2(1,1:end-1));
         enter_var=find(min_value==As2(1,:),1);
         leave_var=0;
         ratio=[];
         for k = 2:height(As2)%%calculate ratio
            if As2(k, end)<=0
                 ratio = [ratio,0];%Put 0 if the aswer will be negative
            elseif As2(k, enter_var) <= 0
                ratio = [ratio,0];
            else
                 ratio = [ratio, As2(k, end) / As2(k, enter_var)];
            end
         end
              min_ratio_value=min(ratio);
              if min_ratio_value==0
                  positive_ratios = ratio(ratio > 0);
                  min_ratio_value=min(positive_ratios);
              end 
              max_ratio_value=max(ratio);
              negative_value=min(As2(1,1:end-1))<0; 
          if negative_value==true && max_ratio_value==0 %%if ratio is all null and we still have negative value in the z vector
              disp('Error this problem is unbound');%%Then this proplem is unbound
              negative_value=false; 
          elseif negative_value==true%%if not the problem continue normally
          leave_var = find(ratio == min_ratio_value,1)+1;
          disp(["leave var=",leave_var]);
          disp(["enter var=",enter_var]);
          pivot=As2(leave_var,enter_var);%%%%%
          disp(['pivot=',num2str(pivot)])
          for g=1:length(As2)
          As2(leave_var,g)=[As2(leave_var,g)/pivot];
          end
          vars_char(leave_var,:)=char(varis(1,enter_var));
          pivot_row=leave_var;
          for t=1:height(As2)
              if t~=pivot_row
                 remainder=As2(t,enter_var);
                As2(t,:)=As2(t,:)-(remainder*As2(leave_var,:));
              end
          end
          end
     end
     disp([vars_char]);
     disp('----The answers to Simplex method----');
     disp(['Z=',num2str(As2(1,end))]);
     for i=2:height(As2)
         disp(['x',num2str(i-1),'=',num2str(As2(i,end))]);
     end
     disp(['Type ratio to check how much was ratio , and vars_char to check the what we call basic variables (the ones on the right)' ...
         'and leave_var to checl leave_var , and enter_var for entering variable' ...
         'and As2 for final result and As1 for the table b4 we start the computation' ...
         'if its 2-phase method use Aa , Aa2 ,Aa3 to check instead of As1 and As2']);