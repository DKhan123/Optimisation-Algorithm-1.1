 function [Result,Status]=GEM(T,P,Bulk)
%% Run meemum.exe 
% Generate executable file for meemum input: Name of dat file (BFCFs in this example), no to
% changing bulk composition, Enter T & P, 0 0 to end meemum. 
Pb=P*1e4;
TP = [T,Pb,Bulk]';
TPn=reshape(TP,size(TP,1)*size(TP,2),1);
% for i=1:length(T)
% TPn(2*i-1)=TP(i,1);
% TPn(2*i)=TP(i,2);
% end

executefile = fopen('executefilename.txt','w');
% Change BFCFS to the Name of Your File
fprintf(executefile,[strcat('BFCFs','\n','y','\n')]);
fprintf(executefile,'%d\n',TPn);
fprintf(executefile,'%s\n','0','0');
fclose(executefile);

jsystem('meemum690.exe   <executefilename.txt   >myFile.txt ');

Result=fileread('myFile.txt');
 end

 
%  
%  tic
%  [a,b]=jsystem('meemum691.exe   <executefilename.txt')
%  toc
