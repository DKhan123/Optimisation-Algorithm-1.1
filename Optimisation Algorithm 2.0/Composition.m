function [Comp_Target,ErrorSd]=Composition(Bulk,Experiments,StablePhases,PhaseData,process)
%%
% Get Bulk Composition
BulkComp=Bulk{1};
BulkCompVal=double(string(BulkComp(:,2:end-1)));
% Get Experiment Composition + Errors
ExpComp=Experiments{1};
% Remove Empty Rows 
ExpComp(sum(ExpComp<=0,2)==length(ExpComp),:)=[];
% Identify Phases Present 
PhaseNames=PhaseData(1,:);
Cations=PhaseData(2,:);
% Wt% or Mol% 
WtMol=string(BulkComp(:,end));
% Identify Components
Component=lower(Bulk{2}(1,:));
Component_Er=lower(Component+"_er");
Components=[Component,Component_Er];
CompData=(ExpComp(:,ismember(lower(Experiments{3}(1:end)),lower(Components))));
% Mr
Mr=str2double(Bulk{2}([2,3],:));

%% Start Calculations 
% n = number of components
n=length(Bulk{2});
% Check for Data 
if ~isempty(CompData)
% Test for Mass Balance     
    if strcmpi(WtMol,'Wt')
        ExpMC=ExpMassCorrection(CompData,BulkCompVal,n);
    else
       CompData(:,1:n)=Mr(1,:).*(CompData(:,1:n)./Mr(2,:));
       ExpMC=ExpMassCorrection(CompData,BulkCompVal,n);
    end
% Convert to Mol                
ExpResults=Mr(2,:).*(ExpMC(:,1:n)./Mr(1,:));
% Cations of Phases 
for i=1:length(StablePhases)
Cats(i)=double(Cations(ismember(PhaseNames,StablePhases(i))));
end

% Normalise Mol to Phase Formula 
ExpMolar=ExpResults./sum(ExpResults,2).*Cats';
ExpError=(CompData(:,n+1:end)./Mr(1,:));
ExpMolarError=ExpError.*Cats';
% Normalise Mol to 1 
ExpTarget=ExpMolar;
ExpTarget(ExpTarget<0)=0;
ExpTarget=ExpTarget./sum(ExpTarget,2).*round(sum(ExpTarget,2));
NormExpMol=ExpTarget./sum(ExpTarget,2);
% Normalise Bulk Comp 
MBulkComp=(BulkCompVal./Mr(1,:).*Mr(2,:))./sum(BulkCompVal./Mr(1,:).*Mr(2,:));
% Reshape Matrix 
l=sum(ExpTarget>-1);
NormExpMol=reshape(NormExpMol(~isnan(NormExpMol)),l(1),n);

% Calculate Mass Fraction
ExpWT=lsqnonneg(NormExpMol',MBulkComp');
ExpWT=ExpWT./sum(ExpWT);
Exp_Target=[ExpWT*100,reshape(ExpTarget(~isnan(ExpTarget)),l(1),n)]';
% Calculate Error 
[~,CI]=regress(MBulkComp',NormExpMol');
WtError=abs(CI(:,1)-CI(:,2))/3.92; % +-1.96 is 95% confidence interval in std 

%% For Stx11 example - Remove Components not present in Thermo Model - Can be Deleted 
if ~ismember(lower(process),'n')
% O polymorph correction 
for i=1:length(StablePhases)
Var=logical(ismember(PhaseNames,StablePhases(i)));
in=find(Var==1);
if any(in==[1,2,3])
Exp_Target([3,6,7],i)=0;
Exp_Target(2,i)=1;
Exp_Target(4,i)=2*Exp_Target(4,i)/(Exp_Target(4,i)+Exp_Target(5,i));
Exp_Target(5,i)=2-Exp_Target(4,i);
end
if any(in==4)
    if Exp_Target(2,i)>2
        Exp_Target(2,i)=2;
    end
    Exp_Target(7,i)=0;
end
if any(in==5)
    if Exp_Target(2,i)>2
        Exp_Target(2,i)=2;
    end
end
if any(in==6)
    Exp_Target([3,6,7],i)=0;
    Exp_Target(2,i)=2;
    Exp_Target(4,i)=2*Exp_Target(4,i)/(Exp_Target(4,i)+Exp_Target(5,i));
    Exp_Target(5,i)=2-Exp_Target(4,i);
end
if any(in==9)
    if Exp_Target(2,i)>4
        Exp_Target(2,i)=4;
        Exp_Target(3:end,i)=8*Exp_Target(3:end,i)/sum(Exp_Target(2:end,i));
    end
end
if any(in==11)
Exp_Target([6,7],i)=0;
Exp_Target(2:5,i)=2*Exp_Target(2:5,i)/sum(Exp_Target(2:5,i));
end
if any(in==12)
Exp_Target([2,3,6,7],i)=0;
Exp_Target(4,i)=Exp_Target(4,i)/(Exp_Target(4,i)+Exp_Target(5,i));
Exp_Target(5,i)=1-Exp_Target(4,i);
end
end
end

% Recalculate 
NormExpMol=Exp_Target(2:end,:)./sum(Exp_Target(2:end,:),1);
ExpWT=lsqnonneg(NormExpMol,MBulkComp');
ExpWT=ExpWT./sum(ExpWT);
[~,CI]=regress(MBulkComp',NormExpMol);
WtError=abs(CI(:,1)-CI(:,2))/3.92;

%% End of Removal 
%%

% Check for Empty Rows
if sum(sum(ExpWT,2)==0)>1
    ExpMC=ExpMassCorrection(ExpComp,BulkCompVal,n);
    ExpResults=Mr(2,:).*(ExpMC(:,1:6)./Mr(1,:));
    ExpMolar=ExpResults./sum(ExpResults,2).*Cats';
    NormExpMolar=(ExpMolar./sum(ExpMolar,2))';
    ExpWT=lsqnonneg(NormExpMolar,MBulkComp');
    ExpWT=ExpWT./sum(ExpWT);
    [~,CI]=regress(MBulkComp',NormExpMolar);
    WtError=abs(CI(:,1)-CI(:,2))/3.92; 
end

Exp_Target(1,:)=ExpWT*100;

% AlO1.5 to Al203 and NaO0.5 to Na2O - For Perple_X
Exp_Target([3,7],:)=Exp_Target([3,7],:)/2;

else
    Exp_Target=[];
end

% Collate 
ErrorSd=[WtError*100,ExpMolarError];
Comp_Target=[Exp_Target;ExpComp(:,ismember(Experiments{3}(2:end),Component_Er))'];
end


