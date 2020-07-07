%% Converts Composition from Wt% to Moles
function [Exp_Target,ExpMolarError]=Composition(BulkComp,ExpR)

% Read Data 
StixData=readtable('StixData.dat');
PhaseNames=string(table2cell(StixData(:,2)))';

% State Phases, Number of cations and Mr of Components 
PhaseNam = ({'O','Wad','Ring','Opx','Cpx','C2/c','ca-pv','Aki','Gt_maj','q','Pv','Wus','Sp','CF','Ppv','Pl'});
Phase_N = [3,3,3,4,4,4,2,3,8,1,2,1,12,3,2,5]; 
Mr=[60.0843,101.961276,71.8444,40.3044,56.0774,61.97894;1,2,1,1,1,2];

[Exp_Phases]=find(sum(ExpR,2)>0);
ExpMC=ExpMassCorrection(ExpR);
ExpTemp=zeros(16,6);
for i = 1:length(Exp_Phases)
    ExpTemp(Exp_Phases(i),:)=ExpMC(i,:);
end
ExpMC=ExpTemp(Exp_Phases,:);
ExpResults=Mr(2,:).*(ExpMC(:,1:6)./Mr(1,:));
ExpMolar=ExpResults./sum(ExpResults,2).*Phase_N(Exp_Phases)';
ExpError=(ExpR(Exp_Phases,7:12)./Mr(1,:));
ExpMolarError=ExpError.*Phase_N(Exp_Phases)';
Exp_Stable=PhaseNam(Exp_Phases);


%% Determine endmember of experiments 
ExpTarget=ExpMolar;
ExpTarget(ExpTarget<0)=0;
ExpTarget=ExpTarget./sum(ExpTarget,2).*round(sum(ExpTarget,2));
NormExpTarget=ExpTarget./sum(ExpTarget,2);
MBulkComp=(BulkComp./Mr(1,:).*Mr(2,:))./sum(BulkComp./Mr(1,:).*Mr(2,:));
l=sum(ExpTarget>-1);
NormExpTarget=reshape(NormExpTarget(~isnan(NormExpTarget)),l(1),6);
ExpWT=lsqnonneg(NormExpTarget',MBulkComp');
ExpWT=ExpWT./sum(ExpWT);
Exp_Target=[ExpWT*100,reshape(ExpTarget(~isnan(ExpTarget)),l(1),6)]';

% Experiment Composition Correction 
for i=1:length(Exp_Stable)
Var=logical(ismember(PhaseNames,Exp_Stable(i)));
in=find(Var==1);
if any(in==1)
Exp_Target([3,6,7],i)=0;
Exp_Target(2,i)=1;
Exp_Target(4,i)=2*Exp_Target(4,i)/(Exp_Target(4,i)+Exp_Target(5,i));
Exp_Target(5,i)=2-Exp_Target(4,i);
end
if any(in==3)
Exp_Target(2,i)=1;
Exp_Target([3,6,7],i)=0;
Exp_Target(4,i)=2*Exp_Target(4,i)/(Exp_Target(4,i)+Exp_Target(5,i));
Exp_Target(5,i)=2-Exp_Target(4,i);
end
if any(in==5)
Exp_Target(2,i)=1;
Exp_Target([3,6,7],i)=0;
Exp_Target(4,i)=2*Exp_Target(4,i)/(Exp_Target(4,i)+Exp_Target(5,i));
Exp_Target(5,i)=2-Exp_Target(4,i);
end
if any(in==7)
    if Exp_Target(2,i)>2
        Exp_Target(2,i)=2;
    end
    Exp_Target(7,i)=0;
end
if any(in==11)
    if Exp_Target(2,i)>2
        Exp_Target(2,i)=2;
    end
end
if any(in==16)
    Exp_Target([3,6,7],i)=0;
    Exp_Target(2,i)=2;
    Exp_Target(4,i)=2*Exp_Target(4,i)/(Exp_Target(4,i)+Exp_Target(5,i));
    Exp_Target(5,i)=2-Exp_Target(4,i);
end
if any(in==22)
    if Exp_Target(2,i)>4
        Exp_Target(2,i)=4;
    end
end
if any(in==30)
Exp_Target([6,7],i)=0;
end
if any(in==33)
Exp_Target([2,3,6,7],i)=0;
Exp_Target(4,i)=Exp_Target(4,i)/(Exp_Target(4,i)+Exp_Target(5,i));
Exp_Target(5,i)=1-Exp_Target(4,i);
end
end

NormExpTarget=Exp_Target(2:end,:)./sum(Exp_Target(2:end,:),1);
l=sum(Exp_Target>-1);
ExpWT=lsqnonneg(NormExpTarget,MBulkComp');
ExpWT=ExpWT./sum(ExpWT);

if sum(sum(ExpWT,2)==0)>1
    ExpMC=ExpMassCorrection(ExpR);
    ExpTemp=zeros(16,6);
for i = 1:length(Exp_Phases)
    ExpTemp(Exp_Phases(i),:)=ExpMC(i,:);
end
    ExpMC=ExpTemp(Exp_Phases,:);
    ExpResults=Mr(2,:).*(ExpMC(:,1:6)./Mr(1,:));
    ExpMolar=ExpResults./sum(ExpResults,2).*Phase_N(Exp_Phases)';
    NormExpMolar=(ExpMolar./sum(ExpMolar,2))';
    ExpWT=lsqnonneg(NormExpMolar,MBulkComp');
    ExpWT=ExpWT./sum(ExpWT);
end

Exp_Target(1,:)=ExpWT*100;

% AlO1.5 to Al203 and NaO0.5 to Na2O
Exp_Target([3,7],:)=Exp_Target([3,7],:)/2;
end 
