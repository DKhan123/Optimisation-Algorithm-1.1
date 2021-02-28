% Experiment Optimum Composition Calculator
function [CorrectedMass,MassBalanceTest] = ExpMassCorr(Exp1,BulkComp,n) 

NMassBulkComp=BulkComp./sum(BulkComp')'*100;

Exp1(sum(Exp1')'==0,:)=[];
NMassExp1=Exp1(:,1:n)./sum(Exp1(:,1:n)')'*100;

Exp1WT=lsqnonneg(NMassExp1',NMassBulkComp');
Error=sum(Exp1WT.*NMassExp1)-NMassBulkComp;

Cor=sum(abs(Exp1WT.*Exp1(:,n+1:end)));
Corr=Error./Cor;

CorrectedMass=NMassExp1-(Corr.*Exp1(:,n+1:end));
CorrectedMass=100*CorrectedMass./sum(CorrectedMass')';
CorrectedMass(CorrectedMass<0)=0;

%% Check
A=lsqnonneg(CorrectedMass',NMassBulkComp');
B=A.*CorrectedMass;
MassBalanceTest=sum((sum(B)-NMassBulkComp))<0.01 && sum((sum(B)-NMassBulkComp))>-0.01;
end
