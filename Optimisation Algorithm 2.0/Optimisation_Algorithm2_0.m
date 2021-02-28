%% Aim: Create Variable list - Variables, Other Variables, Variable Matrix
clc
clear 
format long
%%
%%
% GEM Callsign Determination. Formatted for Perple_X
% Change to the Name of the Thermodynamic Model
Thermomodel='stx11ver.dat';
PerpParam=CallSign(Thermomodel);
%%
%%
% Load Bulk Data 
BulkCompFiles=struct2cell(dir('BulkComp*.dat'));
cellfun(@(x)length(x),BulkCompFiles(1,:));
[~,dex]=sort(cellfun(@(x)length(x),BulkCompFiles(1,:)));
BulkCompFiles=BulkCompFiles(1,dex);
BulkCompData=cell(1,length(BulkCompFiles));
Components=cell(1,length(BulkCompFiles));
BulkId=strings(1,length(BulkCompFiles));
for i=1:length(BulkCompFiles)
    BulkCompData{1,i}=fileread(BulkCompFiles{1,i});
    Components{1,i}=strip(string(split(extractBetween(BulkCompData{1,i},'Components','end'),',')));
    Components{1,i}=reshape(strip(Components{1,i}(1:end-1)),length(Components{1,i}(1:end-1))/3,3)';
    BulkComp(i,:)=cellstr(strip(split(extractBetween(BulkCompData{1,i},'Bulk Comp','end'),',')));
    BulkId(1,i)=BulkComp{i};
end
% Load Experimental Data
DataFiles=struct2cell(dir('Data*.dat'));
cellfun(@(x)length(x),DataFiles(1,:));
[~,dex]=sort(cellfun(@(x)length(x),DataFiles(1,:)));
DataFiles=DataFiles(1,dex);
% Load Phases in Database
format long
Data=readtable('StixData.dat');
[PhaseNames,PId]=unique(Data(:,2),'stable');
PhaseData= [[string(table2cell(PhaseNames)'),"System"];[table2array(Data(PId,32))',0]];
%%
%%
% Check Data Files, Extract P,T Conditions. 
ExpResults=cell(3,length(DataFiles));
T=zeros(1,length(DataFiles));
P=zeros(1,length(DataFiles));
StablePhases=cell(1,length(DataFiles));
for i=1:length(DataFiles)
ExpData=readtable(DataFiles{1,i},'ReadRowNames',true,'HeaderLines',7);
ExpResults{1,i}=table2array(ExpData);
ExpResults(2,i)=ExpData.Properties.DimensionNames(1);
if ExpResults(2,i)=="Row"
    ExpResults{2,i}=BulkId(1);
end
ExpResults{3,i}=string(ExpData.Properties.VariableNames);
NamesStr=string(ExpData.Properties.RowNames)';
WtIn=ismember(lower(ExpResults{3,i}(1,:)),"wt");
WtPc{i}=(ExpResults{1,i}(:,WtIn))';
PTCondition=split(string(extractBetween(fileread(DataFiles{1,i}),'PT Conditions:','End')),',');
PTCondition=reshape(PTCondition(1:4),2,2);
[Namdex,Namein]=ismember(PhaseData(1,:),NamesStr);
StablePhases{i}=PhaseData(1,Namdex);
if sum(contains(PTCondition(1,:),'C','IgnoreCase',true)==1)
T(i)= double(PTCondition(2,contains(PTCondition(1,:),'C','IgnoreCase',true)))+273.15;
else
    T(i)= double(PTCondition(2,contains(PTCondition(1,:),'K','IgnoreCase',true)));
end
if sum(contains(PTCondition(1,:),'GPa','IgnoreCase',true)==1)
P(i)= double(PTCondition(2,contains(PTCondition(1,:),'GPa','IgnoreCase',true)));
else
    P(i)= double(PTCondition(2,contains(PTCondition(1,:),'bar','IgnoreCase',true)));
end
NamesStr=strrep(NamesStr,'_1','');
if sum(Namdex)~= length(unique(NamesStr)) 
    Unknown = [NamesStr(~ismember(NamesStr,StablePhases{i})),i];
    if lower(Unknown{1})==lower('System')
    else
    warning('Unknown Phases Present')
    disp(Unknown)
    break
    end
end
end
StablePhases{i}=NamesStr;
%%
%%
% Process Compositional Data 
Compositions=cell(1,length(T));
Comp_Error=cell(1,length(T));
process=input('Do you wish to process compositional data? Y/N [Y]:','s');
if isempty(process)
    process='Y';
end
Bulk=cell(1,2);
% Format Experimental Data
for i=1:length(T)
   [~,BulkFind]=ismember(ExpResults{2,i},BulkId);
   if BulkFind==0
       BulkFind=1;
   end
Bulk{1}=BulkComp(BulkFind,:);
Bulk{2}=Components{BulkFind,:};
[Compositions{i},Comp_Error{i}]=Composition(Bulk,ExpResults(:,i),StablePhases{i},PhaseData,process);
if sum(WtPc{i})>0
    Compositions{i}(1,:)=[];
    ExpId=ExpResults{3,i}(1,~contains(ExpResults{3,i},'_er'));
elseif sum(ismember(lower(ExpResults{3,i}),lower(Components{1}(1,:))))>1
    ExpId=['Wt',ExpResults{3,i}(1,~contains(ExpResults{3,i},'_er'))];
else 
    ExpId=ExpResults{3,i};
end
 ExpResults{4,i}=Comp_Error{i}(:,2:end)';
 ExpResults{4,i}(ExpResults{4,i}==0)=0.1;
 ExpResults{1,i}=[Compositions{i}(1:end,:)',ExpResults{1,i}(:,~ismember(lower(ExpResults{3,i}(1:end)),lower([Components{1}(1,:),Components{1}(1,:)+'_er'])))]';
 ExpResults{3,i}=ExpId';
 WtError{i}=Comp_Error{i}(:,1)';
end
%%
%%
% Select Parameters To Optimise 
ParamName=string(Data.Properties.VariableNames);
ParamEr=strfind(ParamName,'Er');
ParamIdx = cellfun(@(x)~isempty(x) && x >= 0, ParamEr);
ParamList = [strrep(ParamName(ParamIdx),'_Er',''),'T','P','Wij'];
[ParamOpt,~] = listdlg('ListString',ParamList);
if sum(ismember(ParamOpt,length(ParamOpt)-2:length(ParamOpt)))>1 && sum(ismember(ParamOpt,1:length(ParamOpt)-2))>0
if length(ExpResults)<size(Data,1)
warning('System Underdetermined.')   
end
end
%%
%%
% Format Parameter Values
ParamOptIda=ismember(ParamName,ParamList(ParamOpt));
ParamOptIdb=ismember(ParamName,strcat(ParamList(ParamOpt),'_Er'));
ParamOptId=ParamOptIda+ParamOptIdb;
ParamOptId(1)=[];
DataVal=double(string(table2cell(Data)));
ParamValues=DataVal(:,logical(ParamOptIda));
ParamError=DataVal(:,logical(ParamOptIdb));
%%
%%
% Input Interaction Parameters 
WijFile='WijData.dat';
WijData=readtable(WijFile);
WijName=string(table2array(WijData(:,1)));
Wij=table2array(WijData(:,2));
Wij_Er=table2array(WijData(:,3));
% Set Errors for P,T   
T_Er=ones(1,length(T))*35;
P_Er=ones(1,length(P))*0.1;
%%
%% 
%Collate Parameters
Parameters=[ParamName,'T','P','Wij','T_Er','P_Er','Wij_Er'];
AllVar=[cellstr(Parameters);[num2cell(DataVal,1),T',P',Wij,T_Er',P_Er',Wij_Er]];
ParamYes=ParamList(ParamOpt);
ParamNo=ParamList(~ismember(ParamList,ParamYes));
ErrorYes=strcat(ParamYes,'_Er');
ErrorNo=strcat(ParamNo,'_Er');
Variables=AllVar(2,ismember(AllVar(1,:),ParamYes));
OtherVariables=AllVar(2,ismember(AllVar(1,:),ParamNo));
Comp_Err=AllVar(2,ismember(AllVar(1,:),ErrorYes));
OtherErrors=AllVar(2,ismember(AllVar(1,:),ErrorNo));
VarData=cell2mat(Variables(:));
OtherVarData=cell2mat(OtherVariables(:));
ErrorData=cell2mat(Comp_Err(:));
OtherErrorData=cell2mat(OtherErrors(:));
% Determine Bounds 
UBound=input('Upper bound for optimisation:');
if isempty(UBound)
    UBound=[];
end
LBound=input('Lower bound for optimisation:');
if isempty(LBound)
    LBound=[];
end
LBP=VarData-LBound*ErrorData;
UBP=VarData+UBound*ErrorData;

if any(VarData>0 & LBP<0)
    LBP(VarData>0 & LBP<0)=0;
end

if any(VarData<0 & UBP>0)
    UBP(VarData<0 & UBP>0)=0;
end

%%
%% 
% Calculate the number of each variable
VarSize{1}=[ParamYes;num2cell(cellfun ('length',Variables))];
VarSize{2}=[ParamNo;num2cell(cellfun ('length',OtherVariables))];
VarSize{3}=lower(Components{1}(1,:));
%%

for i =1:length(ExpResults)
    if any(ExpResults{1,i}(ismember(lower(ExpResults{3,i}),"wt"),:)==0)
    if WtError{i}(:,ExpResults{1,i}(ismember(lower(ExpResults{3,i}),"wt"),:)==0)==0
    StablePhases{i}(:,ExpResults{1,i}(ismember(lower(ExpResults{3,i}),"wt"),:)==0)=[];
    ExpResults{4,i}(:,ExpResults{1,i}(ismember(lower(ExpResults{3,i}),"wt"),:)==0)=[];
    WtError{i}(:,ExpResults{1,i}(ismember(lower(ExpResults{3,i}),"wt"),:)==0)=[];
    ExpResults{1,i}(:,ExpResults{1,i}(ismember(lower(ExpResults{3,i}),"wt"),:)==0)=[];
    end
    end
end
%%
%% Compile Data
% Variables, Exponents, Other Variables, Lower Bound, Upper Bound....
% VarInputs=[NormVar,Var10Exp,LBP./(10.^Var10Exp),UBP./(10.^Var10Exp)];
VarInputs=[zeros(length(VarData),1),zeros(length(VarData),1),-LBound*ones(length(VarData),1),UBound*ones(length(VarData),1)];
VarInputs=[VarInputs,VarData,ErrorData];

OtherVariables;
VarSize;
ExperimentData=[ExpResults;StablePhases;num2cell(T);num2cell(P);WtError];
Wij=[WijName,Wij,Wij_Er];
GemParam=[PerpParam;table2cell(Data(:,ismember(Data.Properties.VariableNames,PerpParam(1,:))))];
GemParam(3:end,2:3)=-str2double(GemParam(3:end,2:3));
Compounds=Components{1}(1,:);
%%
StdVar=[];
for i=1:length(Comp_Err)
StdVar=[StdVar;Comp_Err{i}];
end

yn=input('Add covariance data for prior? y/n ','s');
if lower(yn)=='y'
CoV=readtable('Covar.dat');
V0=double(string(table2cell(Cov)));
else
V0=diag(StdVar)^2;
yn='n';
end
%%
i=0;
fv1=0;
while i < 10
[ ~,  ~, fval, x]=VarOutput(VarInputs,VarSize,ExperimentData,Wij,GemParam,BulkComp,0);
VarInputs(:,1)=x;
i=i+1;
if fval==fv1
    break 
end
fv1=fval;
end

[Cov, CovSPD, yminEst, x]=VarOutput(VarInputs,VarSize,ExperimentData,Wij,GemParam,BulkComp,1);

if length(ExpResults)<length(VarData)
else
    yminEst=yminEst/(length(ExpResults)-length(VarData));
end

if yminEst<1
    yminEst=1;
end


StDevParam=sqrt(diag(2*yminEst*Cov));
V1=2*yminEst*Cov;
if sum(real(StDevParam)>0)<length(x)
StDevParam=sqrt(diag(2*yminEst*CovSPD));
V1=2*yminEst*CovSPD;
end

if lower(yn)~='y'
    V1=diag(diag(V1));
end

OptMin=(x.*StdVar)+VarData;
StDevParam=StDevParam.*StdVar;

BayesCov=V1*(inv(V0+V1))*V0;
BayesMean=V1*(inv(V0+V1))*VarData   +  V0*(inv(V0+V1))*OptMin;
BayesStd=sqrt(diag(BayesCov)).*StdVar;

CorrOpt=zeros(length(x));
for i = 1:length(x)
    for j = 1:length(x)
CorrOpt(i,j)=CovSPD(i,j)./(sqrt(CovSPD(i,i))*sqrt(CovSPD(j,j)));
    end
end

Pi=CDF(inf(length(BayesMean),1),BayesMean,BayesCov);
if Pi~=1
    BayesCov=BayesCov./Pi;
    BayesMean=BayesMean./Pi;
end

BC=ones(23,6);
BC=BC.*double(string({'43.68','3.13','18.71','31.50','2.49','0.5'}));

PrintGraph(VarData,OptMin,BayesMean,VarSize,GemParam,Wij,T,P,Data(:,[1,2]),BC)
