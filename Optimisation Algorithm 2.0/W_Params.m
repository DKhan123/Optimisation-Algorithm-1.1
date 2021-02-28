function [DeltaF,W] = W_Param(W,T,P,StablePhases,E_Target,E_Error)
for i =1:length(E_Target)
%1.SD
E_Errora{i}=E_Error{i};
end

ExcessNames=["fo fa)";"wad fwad)";"ring fring)";"odi ts)";"odi en)";"cen di)";"jd  di)";"cts di)";"cen cts)";"cen hed)";"jd  cts)";"hed cfs)";"cfs di)";"aki cor)";"gr maj)";"gr py)";"py maj)";"perov aperov)";"per wus)";"sp herc)";"ppv appv)";"an ab)"];
ExcessW=[7.6;16.5;9.1;48;32.1;24.7;24.3;26;60.6;24.7;10;24.7;24.7;66;58;30;21.3;116;13;5;60;26];
New_W=([ExcessNames,ExcessW]);
SolutionModel='Solution_Model.dat';
FileData=fileread(SolutionModel);
FileData1=FileData;

New_W(:,2)=W;

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
Result=Meemum_O(T,P);

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
%%
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
E_Target{i}=E_Target{i}(:,SE{i});
Error{i}=E_Errora{i}(SE{i},:)'*1;
end
 
%%
for i =1:length(PName)
if length(unique(SPP{i}))~=length(SPP{i})
[IdNam{i},Id{i}]=unique(SPP{i});
Dupe{i}=1:length(SPP{i});
Dupe{i}(Id{i})=[];
DupeID{i} = SPP{i}(Dupe{i});
Duplicate{i}=Predicted{i}(1:end,Dupe{i});
Predicted{i}=Predicted{i}(1:end,Id{i});
SPP{i}(Dupe{i})=[];
SDUP{i}=SPP{i}(:,Dupe{i});
Len_Predict{i}=sum(2*(ismember(SPP{i},StablePhases{i})));

elseif length(unique(SPP{i}))<length(SPP{i})-1
Len_Predict{i}=-2;
else
Len_Predict{i}=sum(2*(ismember(SPP{i},StablePhases{i}))-1);
end
CoPhase=1*sum(ismember(string(SPP{i}),string(SEP{i})));
if ne(Len_Predict{i},length(StablePhases{i}))
    DelTargetA{i} = ((1/CoPhase)*sumabs(Predicted{i}(1,ismember(string(SPP{i}),string(SEP{i})))-E_Target{i}(1,ismember(string(SEP{i}),string(SPP{i})))));
    DelTargetB1{i} = (abs(Predicted{i}(2:end,ismember(string(SPP{i}),string(SEP{i})))-E_Target{i}(2:end,ismember(string(SEP{i}),string(SPP{i})))))-Error{i}(:,ismember(string(SEP{i}),string(SPP{i})));
    DelTargetB{i} = sumabs(max(DelTargetB1{i},0))./sum(sum(Predicted{i}(2:end,ismember(string(SPP{i}),string(SEP{i})))))*100;
    
    if length(unique(SPP{i}))~=length(SPP{i})
        DelDupe{i}=(1*sumabs(Duplicate{i}(1,:))) + sumabs(Duplicate{i}(2:end,:))-E_Target{i}(2:end,Dupe{i}-1)-Error{i}(:,Dupe{i}-1);
    else 
        DelDupe{i}=0;
    end
    
    DelExtra{i} = (1*sumabs(Predicted{i}(1,~ismember(string(SPP{i}),string(SEP{i}))))) + sumabs(Predicted{i}(2:end,~ismember(string(SPP{i}),string(SEP{i}))));
    DelMissing{i} = (1*sumabs(E_Target{i}(1,~ismember(string(SEP{i}),string(SPP{i}))))) + sumabs(E_Target{i}(2:end,~ismember(string(SEP{i}),string(SPP{i}))));
    DelEnd{i}=DelTargetB{i}+DelTargetA{i}+DelExtra{i}+DelMissing{i}+DelDupe{i};
elseif Len_Predict{i}==-2
    DelEnd{i}=100;
else
    DelTargetA{i} = ((1/CoPhase)*sumabs(Predicted{i}(1,ismember(string(SPP{i}),string(SEP{i})))-E_Target{i}(1,ismember(string(SEP{i}),string(SPP{i})))));
    DelTargetB1{i} = (abs(Predicted{i}(2:end,ismember(string(SPP{i}),string(SEP{i})))-E_Target{i}(2:end,ismember(string(SEP{i}),string(SPP{i})))))-Error{i}(:,ismember(string(SEP{i}),string(SPP{i})));
    DelTargetB{i} = sumabs(max(DelTargetB1{i},0))./sum(sum(Predicted{i}(2:end,ismember(string(SPP{i}),string(SEP{i})))))*100;
    DiffEnd{i} = DelTargetB{i}+DelTargetA{i};
    DelEnd{i}=DiffEnd{i};
end

end
for i = 1:length(DelEnd)
    Delta(i)=sum(DelEnd{i}.^2);
end
DeltaF=sum(Delta);
end




