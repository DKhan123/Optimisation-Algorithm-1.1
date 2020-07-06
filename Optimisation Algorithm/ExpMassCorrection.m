% Experiment Optimum Composition Calculator
function [CorrectedMass,MassBalanceTest] = ExpMassCorr(Exp1) 

BulkComp=csvread('BulkComp.dat',1,1);

NMassBulkComp=BulkComp./sum(BulkComp')'*100;

Exp1(sum(Exp1')'==0,:)=[];
NMassExp1=Exp1(:,1:6)./sum(Exp1(:,1:6)')'*100;

Exp1WT=lsqnonneg(NMassExp1',NMassBulkComp');
Error=sum(Exp1WT.*NMassExp1)-NMassBulkComp;

Cor=sum(abs(Exp1WT.*Exp1(:,7:end)));
Corr=Error./Cor;

CorrectedMass=NMassExp1-(Corr.*Exp1(:,7:end));
CorrectedMass=100*CorrectedMass./sum(CorrectedMass')';
CorrectedMass(CorrectedMass<0)=0;

%% Check
A=lsqnonneg(CorrectedMass',NMassBulkComp');
B=A.*CorrectedMass;
MassBalanceTest=sum((sum(B)-NMassBulkComp))<0.001 && sum((sum(B)-NMassBulkComp))>-0.001;
end
