 function [Result,Status]=Meemum_O(T,P)
%% Run meemum.exe 
% Generate executable file for meemum input: Name of dat file (BFCFs in this example), no to
% changing bulk composition, Enter T & P, 0 0 to end meemum. 
Pb=P*1e4;
TP = [T,Pb];
TPn=zeros(numel(TP),1);
for i=1:length(T)
TPn(2*i-1)=TP(i,1);
TPn(2*i)=TP(i,2);
end

executefile = fopen('executefilename.txt','w');
fprintf(executefile,[strcat('BFCFs','\n','n','\n')]);
fprintf(executefile,'%d\n',TPn);
fprintf(executefile,'%s\n','0','0');
fclose(executefile);
[status,Result] = system('meemum.exe<executefilename.txt');
end