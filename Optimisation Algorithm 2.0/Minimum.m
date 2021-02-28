function [DeltaF,Variables] = Minimum2(Var,VarInputs,VarExponents,VarSize,ExpData,Wij,GemParam,BulkComp)
% UnStandardise Variables
Var=[Var;VarInputs(length(Var)+1:end,1)];
Variables=(Var.*VarInputs(:,6))+VarInputs(:,5);

% Assign Variables to GEM Input File 
[T,P]=UpdateModel(VarSize, Variables, GemParam, Wij);

for i = 1:size(BulkComp,1)
    BID(i,1)=string(BulkComp(i,1));
end

for i = 1:length(ExpData)   
    for j = 2:size(BulkComp,2)-1
Bulk(i,j-1)=str2double(cellstr(BulkComp(ismember(BID,ExpData{2,i}),(j))));
    end
end

if isempty(T)
    T=cell2mat(ExpData(6,:)');
end
if isempty(P)
    P=cell2mat(ExpData(7,:)');
end

% Run GEM
Result=GEM(T,P,Bulk);
[GemData,TPCheck]=ReadGEM(Result);

% Check for Numerical Instability 
MissCheck=size(cell2mat(ExpData(6, ~ismember(cell2mat(ExpData(6,:)),TPCheck(:,1)))),2);
% If Unstable Adjust T and Recalulate 
if MissCheck>0&&MissCheck<length(ExpData)
ExpData(6,~ismember(cell2mat(ExpData(6,:)),TPCheck(:,1)))=num2cell(cell2mat(ExpData(6, ~ismember(cell2mat(ExpData(6,:)),TPCheck(:,1))))-0.15);
Result=GEM(cell2mat(ExpData(6,:)'),cell2mat(ExpData(7,:)'),Bulk);
[GemData,~]=ReadGEM(Result);
end

% If GEM Data Cannot be Calculated 
if length(GemData)<length(ExpData)
    DeltaF=NaN;
    return
    
%%
% If GEM Data can be Calculated
else
    % Sort and Format GEM and EXP Data 
for i=1:length(GemData)
    GemData{1,i}=[GemData{1,i};zeros(1,length(GemData{1,i}))];
    ExpData{5,i}=lower(ExpData{5,i}');
    [ExpData{5,i},SortExp]=sort(lower(ExpData{5,i}));
    EXPVAL{i,1}=ExpData{1,i}(:,SortExp)'; %#ok
    ExpSigma{i}=[ExpData{8,i}(:,SortExp);ExpData{4,i}(:,SortExp)]';  %#ok
    GEMVAL{i,1}=[GemData{1,i},GemData{2,i},GemData{3,i}];  %#ok
    GEMVAL{i}(2:end,:)=sortrows(lower(GEMVAL{i}(2:end,:)));  %#ok
    
    % Identify Properties to Optimise for
    OPT{i}=lower(ExpData{3,i}');   %#ok
    OPT{i}(ismember(OPT{i},VarSize{3}))=[];   %#ok
end

% Property ID from GEM
PropIds(1,:)=lower([GemData{1,1}(1,:),GemData{2,1}(1,:),GemData{3,1}(1,:)]);

%%
% Identify Stable Phases 
for i = 1:length(ExpData)
RowId(i,:)=contains(lower(PropIds),lower(ExpData{3,i}'));  %#ok
EPhaseP(1:length(ExpData{5,i}),i)=ExpData{5,i};  %#ok
EPP{i}=ExpData{5,i};  %#ok
GPhaseP(1:length(GEMVAL{i}(2:end,1)),i)=GEMVAL{i}(2:end,1);  %#ok
GPP{i}=GEMVAL{i}(2:end,1);  %#ok
Fields(i)=join(ExpData{5,i});  %#ok

% Identify Possible Redundent Phases Based on Exp Wt% Estimate
if length(ExpData{5,i}(EXPVAL{i,1}(:,1)==0))
Redundent{i}=(ExpData{5,i}(EXPVAL{i,1}(:,1)==0));  %#ok
if any(ismember(GPP{i}(2:end),Redundent{i}))
    Redundent{i}=[];  %#ok
end
else 
    Redundent{i}=[];  %#ok
end

% If Redundent Phase is Not Present in GEM, Remove from Calculations 
% if length(Redundent{i})>0
if ~isempty(Redundent{i})
 EPhaseP((EXPVAL{i,1}(:,1)==0),i)=missing;  %#ok
 EPP{i}(EXPVAL{i,1}(:,1)==0)=[];  %#ok
 ExpData{5,i}(EXPVAL{i,1}(:,1)==0)=[];
 ExpData{1,i}(:,ExpData{1,i}(1,:)==0)=[];
 ExpSigma{i}(EXPVAL{i,1}(:,1)==0,:)=[];
 EXPVAL{i,1}(EXPVAL{i,1}(:,1)==0,:)=[];
end

% Calculate the Difference in No. Phases 
if length(ExpSigma{i})==sum(RowId(i,:))
RowSig(i,:)=contains(lower(PropIds),lower(ExpData{3,i}'));  %#ok
DupePhase(i)=abs(sum(ismember(EPP{i},GPP{i}(2:end)))-sum(ismember(GPP{i}(2:end),EPP{i})));  %#ok
DelPhase1(i)=(DupePhase(i)+length(unique([EPP{i};GPP{i}(2:end)])))-sum(ismember(EPP{i},GPP{i}(2:end)));  %#ok
else
    RowSig(i,:)=RowId(i,:);  %#ok
end
end
DelPhase=DelPhase1';

% List Phases
Phases=rmmissing(unique([EPhaseP;GPhaseP]));
Phases(Phases=='0')='System';

%%
% Format Experimental Data
for i=1:length(ExpData)
for k = 1:length(Phases)
GEMV{k,i}=double(GEMVAL{i}(ismember(GEMVAL{i}(:,1),Phases(k)),:));  %#ok
if any(ismember(ExpData{5,i},Phases(k)))
EXPV{k,i}(1:sum(ismember(ExpData{5,i},Phases(k))),1:length(PropIds))=zeros(sum(ismember(ExpData{5,i},Phases(k))),length(PropIds)); %#ok 
EXPV{k,i}(:,RowId(i,:))=EXPVAL{i}(ismember(ExpData{5,i},Phases(k)),:);  %#ok
SigmaV{k,i}(1:sum(ismember(ExpData{5,i},Phases(k))),1:length(PropIds))=zeros(1,length(PropIds));  %#ok
SigmaV{k,i}(:,RowSig(i,:))=ExpSigma{i}(ismember(ExpData{5,i},Phases(k)),:);  %#ok
else
    EXPV{k,i}(1,1:length(PropIds))=zeros(1,length(PropIds));   %#ok
    EXPV{k,i}(1,:)=NaN(1,length(PropIds));  %#ok
    SigmaV{k,i}(1,1:length(PropIds))=zeros(1,length(PropIds));  %#ok
    SigmaV{k,i}(1,:)=NaN(1,length(PropIds));  %#ok
end
end
end

% Calculate Difference in Properties, Assign S.D. to Correct Exp
A=zeros(length(Phases),length(PropIds),length(ExpData));
B=zeros(length(Phases),length(PropIds),length(ExpData));
D=zeros(length(Phases),length(PropIds),length(ExpData));
for i=1:length(ExpData)
for k = 1:length(Phases)
    % If Duplicate Phase
    if size(EXPV{k,i},1) ~= size(GEMV{k,i},1)&&size(EXPV{k,i},1)*size(GEMV{k,i},1)>0
    A(k,RowId(i,:),i)=sum(abs(EXPV{k,i}(1,RowId(i,:))-GEMV{k,i}(1,RowId(i,:))),1);    
    B(k,RowId(i,:),i)=SigmaV{k,i}(:,RowId(i,:));
    % If Identical Phase
    elseif any(ismember(ExpData{5,i},Phases(k)))&&any(ismember(GEMVAL{i}(:,1),Phases(k)))
    A(k,RowId(i,:),i)=sum(abs(EXPV{k,i}(:,RowId(i,:))-GEMV{k,i}(:,RowId(i,:))),1);
    B(k,RowId(i,:),i)=SigmaV{k,i}(:,RowId(i,:));
    % If Missing Phase
    elseif any(ismember(ExpData{5,i},Phases(k)))&&~any(ismember(GEMVAL{i}(:,1),Phases(k)))
    D(k,RowId(i,:),i)=sum(abs(EXPV{k,i}(1,RowId(i,:))),1); 
    % If Extra Phase
    elseif ~any(ismember(ExpData{5,i},Phases(k)))&&any(ismember(GEMVAL{i}(:,1),Phases(k)))
    D(k,RowId(i,:),i)=sum(abs(GEMV{k,i}(1,RowId(i,:))),1);
    end      
end
end

%%
% Calculate Delta 
% n = S.D. Multiplier (Del/n*Sigma)
n=6;
DelC=zeros(length(ExpData),length(Phases)-1);
% Special Case for Composition 
for k = 2:length(Phases)
for i=1:length(ExpData)
if any(sum(B(k,ismember(lower(GEMVAL{1}(1,:)),VarSize{3}),i))>0)  
    NonZComp=EXPV{k,i}(1,ismember(lower(GEMVAL{1}(1,:)),VarSize{3}))>0;
    DelComp=A(k,ismember(lower(GEMVAL{1}(1,:)),VarSize{3}),i);
    SigmaComp=n*B(k,ismember(lower(GEMVAL{1}(1,:)),VarSize{3}),i);
    DelC(i,k-1)=DelComp(NonZComp)/SigmaComp(NonZComp);

% Set DelC Between 0 and 1. Where if Del=n*Sigma, DelC = 1 
DelC(DelC<0)=0;
DelC(DelC>1)=1;

% Calculate Delta for Other Properties 
for j = 1:length(OPT{i})
DelOpt=ismember(lower(GEMVAL{1}(1,:)),OPT{i}(j));
DelWt{j}(i,k-1)=(A(k,DelOpt,i))./(n*(B(k,DelOpt,i)));
DelWt{j}(DelWt{j}<0)=0;
DelWt{j}(DelWt{j}>1)=1;
end
end
end
end

% Format Results and Collate 
for i  = 1:length(Fields)
    for k = 2:length(Phases)
        for j = 1:length(OPT{i})
        MATCH=EPP{i}(ismember(EPP{i},GPP{i}));
        if ismember(Phases(k),MATCH)
DEL_C(i,k-1)=sum(DelC(i,k-1));
DEL_W{j}(i,k-1)=sum(DelWt{j}(i,k-1));
        else 
            DEL_C(i,k-1)=NaN;
            DEL_W{j}(i,k-1)=NaN;
        end
        end
    end
end

% Calculate Residual for Each Experiment 
DEL_W{j+1}=DEL_C;
for i = 1:length(Fields)
    for j = 1:size(DEL_W,2)
    DEL(i,j)=mean(DEL_W{j}(i,~isnan(DEL_W{j}(i,:))));
    end
end
% Add Non Equal Phase Penalty 
DEL=DEL+DelPhase;

% Sum of Squared Residuals 
Del=DEL.^2;
DeltaF=sum(sum(Del));

if DeltaF<14
    AA=1;
end
if DeltaF<10
    AA=1;
end
if DeltaF<8
    AA=1;
end

end

