function [DeltaF] = Minimise(Variables,OtherVar,Names,StablePhases,E,E_Error,Ti,Pi,Wij)

for i =1:length(E)
%1.SD
E_Errora{i}=E_Error{i};
E_Target{i}=E{i};
end

VNam=[" = ","S0 = ","V0 = ","c1 = ","c2 = ","c3 = ","c4 = ","c5 = ","c6 = ","c7 = ","m0 = ","m1 = "];
list=["F0","K0","Kprime","Dy0","y0","q","S","Sprime","ets0","T","P","Wij"];
[in,dex]=ismember(list,Names{1});
[~,xed]=ismember(list,Names{2});

 if in(1)   
     G0Var=Variables((dex(1)-1)*48+1:(dex(1))*48)*1e3;
 else 
     G0Var=OtherVar(:,3+xed(1))*1e3;
 end
 if in(2)   
     C1Var=Variables((dex(2)-1)*48+1:(dex(2))*48)*1e3;
 else 
     C1Var=OtherVar(:,3+xed(2))*1e3;
 end 
 if in(3)   
     C2Var=Variables((dex(3)-1)*48+1:(dex(3))*48);
 else 
     C2Var=OtherVar(:,3+xed(3));
 end     
 if in(4)   
     C3Var=Variables((dex(4)-1)*48+1:(dex(4))*48);
 else 
     C3Var=OtherVar(:,3+xed(4));
 end 
 if in(5)   
     C4Var=Variables((dex(5)-1)*48+1:(dex(5))*48);
 else 
     C4Var=OtherVar(:,3+xed(5));
 end 
 if in(6)   
     C5Var=Variables((dex(6)-1)*48+1:(dex(6))*48);
 else 
     C5Var=OtherVar(:,3+xed(6));
 end  
 if in(9)   
     C6Var=Variables((dex(9)-1)*48+1:(dex(9))*48);
 else 
     C6Var=OtherVar(:,3+xed(9));
 end   
 if in(7)   
     M0Var=Variables((dex(7)-1)*48+1:(dex(7))*48)*1e4;
 else 
     M0Var=OtherVar(:,3+xed(7))*1e4;
 end   
 if in(8)   
     M1Var=Variables((dex(8)-1)*48+1:(dex(8))*48);
 else 
     M1Var=OtherVar(:,3+xed(8));
 end   
 if in(10)   
     Ti=Variables((dex(10)-1)*48+1:(dex(10)-1)*48+length(E));
 else 
     Ti=Ti';
 end  
 if in(11)   
     Pi=Variables(max(dex(1:9))*48+(dex(10)-1)*48+length(E)+1:(max(dex(1:9))*48+(dex(10)-1)*48+2*length(E)));
 else 
     Pi=Pi';
 end
 if in(12)   
     Wij=Variables((dex(12)-1)*48+1:(dex(12)-1)*48+1+length(Wij));
 else 
     Wij=Wij;
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

NewParam=fopen('Stx11_New.dat','w');
fprintf(NewParam,[strcat(Stx11verN)]);
fclose(NewParam);

if in(12)
ExcessNames=["fo fa)";"wad fwad)";"ring fring)";"odi ts)";"odi en)";"cen di)";"jd  di)";"cts di)";"cen cts)";"cen hed)";"jd  cts)";"hed cfs)";"cfs di)";"aki cor)";"gr maj)";"gr py)";"py maj)";"perov aperov)";"per wus)";"sp herc)";"ppv appv)";"an ab)"];
ExcessW=[7.6;16.5;9.1;48;32.1;24.7;24.3;26;60.6;24.7;10;24.7;24.7;66;58;30;21.3;116;13;5;60;26];
New_W=([ExcessNames,ExcessW]);
SolutionModel='Solution_Model.dat';
FileData=fileread(SolutionModel);
FileData1=FileData;

New_W(:,2)=Wij;

New_WJ=join(New_W);
FileDataString1=string(FileData1);
ExtractW1=extractBetween(FileDataString1,'begin_excess_function','end_excess_function');
WReplaceTxt1=extractBetween(FileDataString1,'W(','d3');
        
for j = 1:size(WReplaceTxt1,1)
    FileData1=strrep(FileData1,WReplaceTxt1(j),New_WJ(j));
end
        
test=fopen('Solution_Model_1.dat','w');
fprintf(test,[strcat(FileData1)]);
fclose(test);
Result=Meemum_O(Ti,Pi);
else
    Result=Meemum_O(Ti,Pi);
end

%% Phase Data Extract
PComp=extractBetween(Result,'Phase Compositions (molar  proportions):','Phase speciation');
PhaseComposition=erase(PComp,'%');

for i=1:length(PComp)
PhaseComp{i}=(split(PhaseComposition{i}));
PhaseNum{i}=single(str2double(PhaseComp{i}));
n(i)=length(PhaseComp{i})/11;
for j=1:n(i)
 Phase_Data{i}(:,j)=num2cell(PhaseNum{i}((j*11)-9:j*11));
 Phase_Id{i}(j)=PhaseComp{i}(11*j-10);
end
Phase_Data{i}(:,1)=PhaseComp{i}([2:11],1);
Phase_Data{i}=vertcat(Phase_Id{i},Phase_Data{i});
end
%% PhaseSpec Extract
PSpec=extractBetween(Result,'Molar Properties and Density:','Seismic Properties:');

for i=1:length(PSpec)
    PhaseSpec{i}=(split(PSpec{i}));
    PhaseSpec{i}(end)=[];
    PhaseSpec{i}=reshape(PhaseSpec{i},10,numel(PhaseSpec{i})/10);
    VSpec{i}=PhaseSpec{i}(5,2:end-1);
end
%Predicted{i}=[];    
%% Del
for i=1:length(Phase_Data)
Predicted{i}=cell2mat(Phase_Data{i}([2,6:end],2:end));
PName{i}=Phase_Data{i}(1,2:end);
[SPP{i},SP{i}]=sort(PName{i});

if any (E_Target{i}(1,:)==0)
Clear{i}=(E_Target{i}(1,:)==0);
E_Target{i}(:,logical(Clear{i}))=[];
StablePhases{i}(logical(Clear{i}))=[];
E_Errora{i}(logical(Clear{i}),:)=[];
end

[SEP{i},SE{i}]=sort(StablePhases{i});
Predicted{i}=Predicted{i}(:,SP{i});
VSpec{i}=VSpec{i}(:,SP{i});
E_Target{i}=E_Target{i}(:,SE{i});
Error{i}=E_Errora{i}(SE{i},:)'*1;
end
 
%%
for i =1:length(PName)
if length(unique(SPP{i}))==length(SPP{i})-1
[~,Id{i}]=unique(SPP{i});
DuplicateId{i}=(sum(1:length(SPP{i}))-sum(Id{i}));
MergeDuplicate{i}=mean([Predicted{i}(2:end,DuplicateId{i}),Predicted{i}(2:end,DuplicateId{i}-1)],2);
MergeWt{i}=Predicted{i}(1,DuplicateId{i})+Predicted{i}(1,DuplicateId{i}-1);
Predicted{i}(:,DuplicateId{i}-1)=[MergeWt{i};MergeDuplicate{i}];
Predicted{i}(:,DuplicateId{i})=[];
%PName{i}(DuplicateId{i})=[];
SPP{i}(DuplicateId{i})=[];
Len_Predict{i}=sum(2*(ismember(string(SPP{i}),string(StablePhases{i})))-1);
elseif length(unique(SPP{i}))<length(SPP{i})-1
Len_Predict{i}=-2;
else
Len_Predict{i}=sum(2*(ismember(SPP{i},StablePhases{i}))-1);
end

if ne(Len_Predict{i},length(StablePhases{i}))
    DelTargetA{i} = (0.25*sumabs(Predicted{i}(1,ismember(string(SPP{i}),string(SEP{i})))-E_Target{i}(1,ismember(string(SEP{i}),string(SPP{i})))));
    DelTargetB{i} = (abs(Predicted{i}(2:end,ismember(string(SPP{i}),string(SEP{i})))-E_Target{i}(2:end,ismember(string(SEP{i}),string(SPP{i})))))-Error{i}(:,ismember(string(SEP{i}),string(SPP{i})));
    DelTargetB{i} = sumabs(max(DelTargetB{i},0));
   % DelTarget{i} = DelTargetA{i}+DelTargetB{i};

    DelExtra{i} = (1*sumabs(Predicted{i}(1,~ismember(string(SPP{i}),string(SEP{i}))))) + sumabs(Predicted{i}(2:end,~ismember(string(SPP{i}),string(SEP{i}))));
    DelMissing{i} = (1*sumabs(E_Target{i}(1,~ismember(string(SEP{i}),string(SPP{i}))))) + sumabs(E_Target{i}(2:end,~ismember(string(SEP{i}),string(SPP{i}))));
    DelEnd{i}=DelTargetB{i}+DelTargetA{i}+DelExtra{i}+DelMissing{i};
elseif Len_Predict{i}==-2
    DelEnd{i}=30;
else
    DelTargetA{i} = (0.25*sumabs(Predicted{i}(1,ismember(string(SPP{i}),string(SEP{i})))-E_Target{i}(1,ismember(string(SEP{i}),string(SPP{i})))));
    DelTargetB{i} = (abs(Predicted{i}(2:end,ismember(string(SPP{i}),string(SEP{i})))-E_Target{i}(2:end,ismember(string(SEP{i}),string(SPP{i})))))-Error{i}(:,ismember(string(SEP{i}),string(SPP{i})));
    DelTargetB{i} = sumabs(max(DelTargetB{i},0));
    DiffEnd{i} = DelTargetB{i}+DelTargetA{i};
    DelEnd{i}=DiffEnd{i};
end

end
for i = 1:length(DelEnd)
    Delta(i)=sum(DelEnd{i}.^2);
end
DeltaF=sum(Delta);
end




