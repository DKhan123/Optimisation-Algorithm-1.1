% Mass Balance Calculator
function [CorrectedMass,MassBalanceTest] = ExpMassCorr(Exp1) 

% Read Bulk Composition File (Wt%)
BulkComp=csvread('BulkComp.dat',1,1);

NMassBulkComp=BulkComp./sum(BulkComp')'*100;

% Remove Excess Phases
Exp1(sum(Exp1')'==0,:)=[];

NMassExp1=Exp1(:,1:6)./sum(Exp1(:,1:6)')'*100;

% Mass Balance (non negative least squares method)
Exp1WT=lsqnonneg(NMassExp1',NMassBulkComp');

% Calculate Error
Error=sum(Exp1WT.*NMassExp1)-NMassBulkComp;

% Correct For Misfit
Cor=sum(abs(Exp1WT.*Exp1(:,7:end)));
Corr=Error./Cor;
CorrectedMass=NMassExp1-(Corr.*Exp1(:,7:end));
CorrectedMass=100*CorrectedMass./sum(CorrectedMass')';
CorrectedMass(CorrectedMass<0)=0;

%% Check
ToleranceCheck=0.001;
A=lsqnonneg(CorrectedMass',NMassBulkComp');
B=A.*CorrectedMass;
MassBalanceTest=sum((sum(B)-NMassBulkComp))<ToleranceCheck && sum((sum(B)-NMassBulkComp))>-ToleranceCheck;
end
