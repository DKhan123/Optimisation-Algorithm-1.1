format long
%% 
% Input Temperature and Pressure for Experiments. 
T=[1100,1200,1250,1250,1400,1450,1500,1550,1575,1650,1700,1710,1725,1750,1750,1750,1750,1750,1765]+273.15;
P=[2,3,4,5,6,7,8,9,10,11,13,13.5,14,15,17,19,20,21,23.5];
% Input Bulk Comp and Data File Names. 
BulkComp=csvread('BulkComp.dat',1,1);
DataFiles=struct2cell(dir('Data*.dat'));
cellfun(@(x)length(x),DataFiles(1,:));
[~,dex]=sort(cellfun(@(x)length(x),DataFiles(1,:)));
DataFiles=DataFiles(1,dex);
for i=1:length(DataFiles)
ExpResults{1,i}=csvread(DataFiles{1,i},1,1);
end

PhaseNam = ({'O','Wad','Ring','Opx','Cpx','C2/c','ca-pv','Aki','Gt_maj','q','Pv','Wus','Sp','CF','Ppv','Pl'});

%%
Exp_Phases=cell(1,length(T));
StablePhases=cell(1,length(T));
E=cell(1,length(T));
E_Error=cell(1,length(T));

for i=1:length(T)
[Exp_Phases{i}]=find(sum(ExpResults{i},2)>0);
StablePhases{i}=PhaseNam(Exp_Phases{i});
[E{i},E_Error{i}]=Composition(BulkComp,ExpResults{i});
end

Data=readtable('StixData.dat');
ParamName=string(Data.Properties.VariableNames);
ParamEr=strfind(ParamName,'Er');
ParamIdx = cellfun(@(x)~isempty(x) && x >= 0, ParamEr);
list = [strrep(ParamName(ParamIdx),'_Er',''),'T','P','Wij'];
[ParamOpt,~] = listdlg('ListString',list);

if sum(ParamOpt>9)>0 && sum(ParamOpt<9)>0
if length(ExpResults)<size(Data,1)
warning('System Underdetermined.')   
CheckParam=input('Solve for T/P/Wij iteratively? Y/N [Y]','s');
if isempty(CheckParam)
   CheckParam = 'Y';
end
end
else
    CheckParam = 'N';
end
MultiStart = 'N';

%% if MultiStart == 'Y'
%% MultiSeed=input('Number of seed values?')
%% end

ParamOptIda=ismember(ParamName,list(ParamOpt));
ParamOptIdb=ismember(ParamName,strcat(list(ParamOpt),'_Er'));
ParamOptId=ParamOptIda+ParamOptIdb;
ParamOptId(1)=[];
Name=string(table2cell(Data(:,1)))';
PhaseNames=string(table2cell(Data(:,2)))';
DataN=double(string(table2cell(Data(:,2:end))));
DataN(:,[5,6,20,21])=[DataN(:,5)/1e3,DataN(:,6)*10,DataN(:,20)/1e3,DataN(:,21)*10];
ParamAll=DataN(:,logical(ParamOptId));
ParamValues=ParamAll(:,1:size(ParamAll,2)/2);
ParamError=ParamAll(:,(size(ParamAll,2)/2)+1:end);
Variables=ParamValues(:);
VarError=ParamError(:);
Variables0=Variables;
VarError0=VarError;

%%
ExcessNames=["fo fa)";"wad fwad)";"ring fring)";"odi ts)";"odi en)";"cen di)";"jd  di)";"cts di)";"cen cts)";"cen hed)";"jd  cts)";"hed cfs)";"cfs di)";"aki cor)";"gr maj)";"gr py)";"py maj)";"perov aperov)";"per wus)";"sp herc)";"ppv appv)";"an ab)"];
T_Error=ones(length(T),1)*35/3;
P_Error=ones(length(T),1)*0.1/3;
W=[7.6;16.5;9.1;48;32.1;24.7;24.3;26;60.6;24.7;10;24.7;24.7;66;58;30;21.3;116;13;5;60;26];
W_Error=[2.2;2.2;2.2;11;1;2;2;4;8.8;2;4;2;2;10;17;5;6.5;10;1;5;10;2];
W0=W;

ParamOptIdc=~ismember(ParamName([6:11,14:16]),list(ParamOpt));
Paramlogic=logical(zeros(1,size(DataN,2)));
Paramlogic([6:11,14:16])=ParamOptIdc;
OtherVarNames=ParamName(Paramlogic);
VariableNames=list(ParamOpt);
Names={VariableNames,OtherVarNames};
Paramlogic(1)=[];
OtherVar=[-DataN(:,11),-DataN(:,4),DataN(:,17),DataN(:,Paramlogic)];

if sum(ParamOpt>9)>0 
    if CheckParam ~= 'Y' %% &&  MultiStart=='Y'
        if sum(ParamOpt==10)>0
            Variables=[Variables;T'];
            VarError=[VarError;T_Error];
        end
        if sum(ParamOpt==11)>0
            Variables=[Variables;P'];
            VarError=[VarError;P_Error];
        end    
        if sum(ParamOpt==12)>0
            Variables=[Variables;W];
            VarError=[VarError;W_Error];
        end
    elseif CheckParam == 'Y'
        ParamOpt(ParamOpt>9)=[];
        VariableNames=list(ParamOpt);
        Names={VariableNames,OtherVarNames};
    end
end

%% if MultiStart=='Y'
%% SeedVal=randn(length(Variables),MultiSeed).*VarError+Variables;
%% end
LBP=Variables-VarError*3;
UBP=Variables+VarError*3;
LBW=W-3*W_Error;
UBW=W+3*W_Error;

SolutionModel='Solution_Model.dat';
FileData=fileread(SolutionModel);
SModel=fopen('Solution_Model_1.dat','w');
fprintf(SModel,[strcat(FileData)]);
fclose(SModel);

Iterate = 0;
while Iterate<5
tic
[NewVariable,F_VAL,History]=VariableOut(Variables,OtherVar,Names,StablePhases,E,E_Error,LBP,UBP,T,P,W);
toc

VNam=[" = ","S0 = ","V0 = ","c1 = ","c2 = ","c3 = ","c4 = ","c5 = ","c6 = ","c7 = ","m0 = ","m1 = "];
list=["F0","K0","Kprime","Dy0","y0","q","S","Sprime","ets0","T","P","Wij"];
[in,dex]=ismember(list,Names{1});
[~,xed]=ismember(list,Names{2});

 if in(1)   
     G0Var=NewVariable((dex(1)-1)*48+1:(dex(1))*48)*1e3;
 else 
     G0Var=OtherVar(:,3+xed(1))*1e3;
 end
 if in(2)   
     C1Var=NewVariable((dex(2)-1)*48+1:(dex(2))*48)*1e3;
 else 
     C1Var=OtherVar(:,3+xed(2))*1e3;
 end 
 if in(3)   
     C2Var=NewVariable((dex(3)-1)*48+1:(dex(3))*48);
 else 
     C2Var=OtherVar(:,3+xed(3));
 end     
 if in(4)   
     C3Var=NewVariable((dex(4)-1)*48+1:(dex(4))*48);
 else 
     C3Var=OtherVar(:,3+xed(4));
 end 
 if in(5)   
     C4Var=NewVariable((dex(5)-1)*48+1:(dex(5))*48);
 else 
     C4Var=OtherVar(:,3+xed(5));
 end 
 if in(6)   
     C5Var=NewVariable((dex(6)-1)*48+1:(dex(6))*48);
 else 
     C5Var=OtherVar(:,3+xed(6));
 end  
 if in(9)   
     C6Var=NewVariable((dex(9)-1)*48+1:(dex(9))*48);
 else 
     C6Var=OtherVar(:,3+xed(9));
 end   
 if in(7)   
     M0Var=NewVariable((dex(7)-1)*48+1:(dex(7))*48)*1e4;
 else 
     M0Var=OtherVar(:,3+xed(7))*1e4;
 end   
 if in(8)   
     M1Var=NewVariable((dex(8)-1)*48+1:(dex(8))*48);
 else 
     M1Var=OtherVar(:,3+xed(8));
 end   
 S0=OtherVar(:,1);
 V0=OtherVar(:,2);
 C7Var=OtherVar(:,3);
Var=[G0Var,S0,V0,C1Var,C2Var,C3Var,C4Var,C5Var,C6Var,C7Var,M0Var,M1Var];

 for i=1:12
     for j=1:48
     VCom(j,i)=join([VNam(i),num2str(Var(j,i))]);
     end
 end
 
for i=1:48
VStruct1(i)=join(VCom(i,1:3));
VStruct2(i)=join(VCom(i,4:10));
VStruct3(i)=join(VCom(i,11:12));
VStruct(i)=VStruct1(i) + newline + VStruct2(i) + newline + VStruct3(i) + newline ;
end

Stx11ver=fileread('stx11ver.dat');

Stx11verA=eraseBetween(Stx11ver,1,"end ");
Stx11verA=eraseBetween(Stx11verA,1,145);

OldParam=extractBetween(Stx11verA,'G0','end');
Stx11verN=Stx11ver;

for i=1:48
Stx11verN=strrep(Stx11verN,OldParam(i),VStruct(i));
end
Variables=NewVariable;
NewParam=fopen('Stx11_New.dat','w');
fprintf(NewParam,strcat(Stx11verN));
fclose(NewParam);

%%
if CheckParam == 'Y'
Iterate=Iterate+1;

tic
[NewW,W_fval,W_history] = Wout(W,T',P',StablePhases,E,E_Error,LBW,UBW);
toc
W=NewW;
New_W=([ExcessNames,W]);
New_WJ=join(New_W);
SolutionModel='Solution_Model.dat';
FileData=fileread(SolutionModel);
FileData1=FileData;
FileDataString1=string(FileData1);
ExtractW1=extractBetween(FileDataString1,'begin_excess_function','end_excess_function');
WReplaceTxt1=extractBetween(FileDataString1,'W(','d3');
        
for j = 1:size(WReplaceTxt1,1)
    FileData1=strrep(FileData1,WReplaceTxt1(j),New_WJ(j));
end
        
SModel=fopen('Solution_Model_1.dat','w');
fprintf(SModel,strcat(FileData1));
fclose(SModel);

%MultiVal(:,MultiSeed)=[mean(History,2),mean(W_history,2)];
else 
    Iterate=5;
end
end

if CheckParam=='Y'
Variables=[Variables0;W0];
NewVariable=[NewVariable;NewW];
VarError=[VarError0;W_Error];
LBP=[LBP;LBW];
UBP=[UBP;UBW];
end

MeanStix=Variables;
StixSdev=VarError;
MeanParam=NewVariable;
if MultiStart=='N'
SdevParam=StixSdev;
else 
SdevParam
end

figure;
for i=1:length(MeanStix)
    [BayesParam(i,1),BayesSdev(i,1)]=Bayes(MeanStix(i),MeanParam(i),SdevParam(i),StixSdev(i));
end
%F0 in KJ

%% Write files
%[OptimProperties,BayesProperties,StixProperties]=WriteFile(NewVariable,NewW,BayesParam,OtherVar,T',P');


