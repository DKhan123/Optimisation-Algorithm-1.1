% Gibbs Free Energy Minimiser 
% This function can be replaced with an alternative Gibbs energy minimiser. 
% If replaced, read files must be updated. 
function [Result,Status]=Meemum_O(T,P)

% Format Pressure and Temperature to relevant units for minimiser [bars, K]
Pb=P*1e4;
TP = [T,Pb];
TPn=zeros(numel(TP),1);
for i=1:length(T)
TPn(2*i-1)=TP(i,1);
TPn(2*i)=TP(i,2);
end

% Generate executable file for meemum input
executefile = fopen('executefilename.txt','w');
% Input commands requested by minimiser: Name of dat file (BFCFs in this
% example), no to changing bulk composition, Enter T & P, 0 0 to end meemum. 
fprintf(executefile,[strcat('BFCFs','\n','n','\n')]);
fprintf(executefile,'%d\n',TPn);
fprintf(executefile,'%s\n','0','0');
fclose(executefile);
[status,Result] = system('meemum.exe<executefilename.txt');
end
