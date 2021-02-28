function[OptimProperties,BayesProperties,StixProperties]=WriteFile(NewVariable,NewW,BayesParam,OtherVar,T,P)

%% Optimised Params
G0Var=NewVariable(1:48)*1e3;
C1Var=NewVariable(49:96)*1e3;
C2Var=NewVariable(97:144);
C3Var=NewVariable(145:192);
C4Var=NewVariable(193:240);
C5Var=NewVariable(241:288);
S0=OtherVar(:,1);
V0=OtherVar(:,2);
C6Var=OtherVar(:,6);
C7Var=OtherVar(:,3);
n=(length(NewVariable)-288)/2;
Ti=NewVariable(289:288+n);
Pi=NewVariable(289+n:end);
NewVar=[G0Var,S0,V0,C1Var,C2Var,C3Var,C4Var,C5Var,C6Var,C7Var];

ExcessNames=["fo fa)";"wad fwad)";"ring fring)";"odi ts)";"odi en)";"cen di)";"jd  di)";"cts di)";"cen cts)";"cen hed)";"jd  cts)";"hed cfs)";"cfs di)";"aki cor)";"gr maj)";"gr py)";"py maj)";"perov aperov)";"per wus)";"sp herc)";"ppv appv)";"an ab)"];
ExcessW=[7.6;16.5;9.1;48;32.1;24.7;24.3;26;60.6;24.7;10;24.7;24.7;66;58;30;21.3;116;13;5;60;26];
New_W=([ExcessNames,ExcessW]);
New_W(:,2)=NewW;

%%
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
        
test=fopen('OptimSolutionModel.dat','w');
fprintf(test,[strcat(FileData1)]);
fclose(test);

VNam=[" = ","S0 = ","V0 = ","c1 = ","c2 = ","c3 = ","c4 = ","c5 = ","c6 = ","c7 = "];
VCom=strings(48,10);

 for i=1:10
     for j=1:48
     VCom(j,i)=join([VNam(i),num2str(NewVar(j,i))]);
     end
 end
 
VStruct1=strings(1,48); 
VStruct2=strings(1,48); 
VStruct=strings(1,48); 

for i=1:48
VStruct1(i)=join(VCom(i,1:3));
VStruct2(i)=join(VCom(i,4:end));
VStruct(i)=VStruct1(i) + newline + VStruct2(i) + newline ;
end

Stx11ver=fileread('Stx11ver.dat');

Stx11verA=eraseBetween(Stx11ver,1,"end ");
Stx11verA=eraseBetween(Stx11verA,1,145);

OldParam=extractBetween(Stx11verA,'G0','m0');
Stx11verN=Stx11ver;

for i=1:48
Stx11verN=strrep(Stx11verN,OldParam(i),VStruct(i));
end

NewParam=fopen('Stx11_OptimParam.dat','w');
fprintf(NewParam,strcat(Stx11verN));
fclose(NewParam);


%% Bayes Params
G0Var=BayesParam(1:48)*1e3;
C1Var=BayesParam(49:96)*1e3;
C2Var=BayesParam(97:144);
C3Var=BayesParam(145:192);
C4Var=BayesParam(193:240);
C5Var=BayesParam(241:288);
S0=OtherVar(:,1);
V0=OtherVar(:,2);
C6Var=OtherVar(:,6);
C7Var=OtherVar(:,3);
BayesW=BayesParam(289:end);
NewVar=[G0Var,S0,V0,C1Var,C2Var,C3Var,C4Var,C5Var,C6Var,C7Var];

New_W=([ExcessNames,ExcessW]);
New_W(:,2)=BayesW;

%%
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
        
test=fopen('BayesSolutionModel.dat','w');
fprintf(test,[strcat(FileData1)]);
fclose(test);

VNam=[" = ","S0 = ","V0 = ","c1 = ","c2 = ","c3 = ","c4 = ","c5 = ","c6 = ","c7 = "];
VCom=strings(48,10);

 for i=1:10
     for j=1:48
     VCom(j,i)=join([VNam(i),num2str(NewVar(j,i))]);
     end
 end
 
VStruct1=strings(1,48); 
VStruct2=strings(1,48); 
VStruct=strings(1,48); 

for i=1:48
VStruct1(i)=join(VCom(i,1:3));
VStruct2(i)=join(VCom(i,4:end));
VStruct(i)=VStruct1(i) + newline + VStruct2(i) + newline ;
end

Stx11ver=fileread('Stx11ver.dat');

Stx11verA=eraseBetween(Stx11ver,1,"end ");
Stx11verA=eraseBetween(Stx11verA,1,145);

OldParam=extractBetween(Stx11verA,'G0','m0');
Stx11verN=Stx11ver;

for i=1:48
Stx11verN=strrep(Stx11verN,OldParam(i),VStruct(i));
end

NewParam=fopen('Stx11_BayesParam.dat','w');
fprintf(NewParam,strcat(Stx11verN));
fclose(NewParam);

%% Results 
StixResult=Meemum_N('BF',T,P);
BayesResult=Meemum_N('BayesBF',T,P);
OptimResult=Meemum_N('OptimBF',T,P);
OB={OptimResult;BayesResult;StixResult};
%% PhaseData
for obn=1:3
PComp=extractBetween(OB{obn},'Phase Compositions (molar  proportions):','Phase speciation');
PhaseComposition=erase(PComp,'%');
Phase_Data=[];
Phase_Id=[];
PhaseComp=cell(1,length(T));
PhaseNum=cell(1,length(T));
Phase_Data=cell(1,length(T));
Phase_Id=cell(1,length(T));

for i=1:length(PComp)
PhaseComp{i}=(split(PhaseComposition{i}));
PhaseNum{i}=single(str2double(PhaseComp{i}));
n(i)=length(PhaseComp{i})/11;
for j=1:n(i)
 Phase_Data{i}(:,j)=num2cell(PhaseNum{i}((j*11)-9:j*11));
 Phase_Id{i}(j)=PhaseComp{i}(11*j-10);
end
Phase_Data{i}(:,1)=PhaseComp{i}(2:11,1);
Phase_Data{i}=vertcat(Phase_Id{i},Phase_Data{i});
if obn ==1
    OptimPhaseData{i}=Phase_Data{i};
elseif obn==2
    BayesPhaseData{i}=Phase_Data{i};
else 
    StixPhaseData{i}=Phase_Data{i};
end
end

%% Speciation Data
Name=["fo","fa","wad","fwad","ring","fring","en","fs","ts","odi","di","hed","cen","cts","jd","c2/c","fc2/c","ca-pv","aki","faki","cor","py","alm","gr","maj","jmaj","q","coe","st","perov","fperov","aperov","per","wus","cfs","sp","herc","an","ab","seif","ppv","fppv","appv","ky","mfer","ffer","nfer","neph"];
PSpec=extractBetween(OB{obn},'Phase speciation (molar proportions):','Molar Properties');
PSpec=erase(PSpec,',');
for i=1:length(PSpec)
PhaseSpeciation{i}=split(PSpec{i});
PhaseSpeciation{i}=PhaseSpeciation{i}(3:end-1);
PhaseSpeciation{i}=erase(PhaseSpeciation{i},':');

PSpecNam{i}=PhaseSpeciation{i}(find(ismember(PhaseSpeciation{i},Name)));
Speciation{i}=single(str2double(PhaseSpeciation{i}(1+find(ismember(PhaseSpeciation{i},Name)))));
Cell_Spec{i}=[PSpecNam{i},num2cell(Speciation{i})];

if obn ==1
    OptimCellSpec{i}=Cell_Spec{i};
elseif obn==2
    BayesCellSpec{i}=Cell_Spec{i};
else
    StixCellSpec{i}=Cell_Spec{i};
end
end

%% Molar Properties
PMol=extractBetween(OB{obn},'Molar Properties and Density:','Seismic Properties:');
for i=1:length(PMol)
    PMolar{i}=split(PMol{i});
    PMolar{i}(end)=[];
    MolProp{i}=reshape(PMolar{i},10,length(PMolar{i})/10)';

if obn ==1
    OptimMolProp{i}=MolProp{i};
elseif obn==2
    BayesMolProp{i}=MolProp{i};
else
    StixMolProp{i}=MolProp{i};
end
end

%% Seismic Properties
PSeis=extractBetween(OB{obn},'Seismic Properties:','Bulk Composition:');
for i= 1:length(PSeis)
    PSeism{i}=split(PSeis{i});
    PSeism{i}(end)=[];
    PSeism{i}(9)=[];
    SeisProp{i}=reshape(PSeism{i},8,length(PSeism{i})/8)';
    
if obn ==1
    OptimSeisProp{i}=SeisProp{i};
elseif obn==2
    BayesSeisProp{i}=SeisProp{i};
else 
    StixSeisProp{i}=SeisProp{i};
end
end

end
Properties=cellstr(["Phase Data";"Speciation";"Molar Properties";"Seismic Properties"]);
OptimProperties=[Properties,[OptimPhaseData;OptimCellSpec;OptimMolProp;OptimSeisProp]];
BayesProperties=[Properties,[BayesPhaseData;BayesCellSpec;BayesMolProp;BayesSeisProp]];
StixProperties=[Properties,[StixPhaseData;StixCellSpec;StixMolProp;StixSeisProp]];

end


