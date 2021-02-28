function [T,P]= UpdateModel(VarSize, Variables, GemParam, Wij)
for i=1:size(VarSize{1},2)
    NewGemParam=Variables(sum(double(VarSize{1,1}(2,1:i-1)))+1:sum(double(VarSize{1,1}(2,1:i))));
    if sum(ismember(GemParam(1,:),VarSize{1,1}(1,i)))>0
    GemParam(3:end,ismember(GemParam(1,:),VarSize{1,1}(1,i)))=(num2str(NewGemParam,9));
    end
    if ismember(VarSize{1}(1,i),'T')
        T=NewGemParam;
    else T=[];
    end
    if ismember(VarSize{1}(1,i),'P')
        P=NewGemParam;
    else
        P=[];
    end
end

VCom=strings(size(GemParam,1)-2,size(GemParam,2));
for j=1:size(GemParam,2)
    for i=1:size(GemParam,1)-2
    VCom(i,j)=join([strip(GemParam(2,j)),GemParam(i+2,j)]);
    end
end

VStruct=string(size(GemParam,1)-2);
for i=1:size(GemParam,1)-2
VStruct1(i)=join(VCom(i,1:3));
VStruct2(i)=join(VCom(i,4:10));
VStruct3(i)=join(VCom(i,11:12));
VStruct(i)=VStruct1(i) + newline + VStruct2(i) + newline + VStruct3(i) + newline ;
end

% Read Original Data File
Stx11ver=fileread('stx11ver.dat');
%Extract Old Parameters
Stx11ver2=eraseBetween(Stx11ver,1,"end ");
Stx11ver2=eraseBetween(Stx11ver2,1,145);
OldParam=extractBetween(Stx11ver2,'G0','end');
% Remove Extra Text
if any(contains(OldParam,'transition'))
  Other = OldParam(contains(OldParam,'transition'));
  OldParam(contains(OldParam,'transition'))=extractBetween(Other,1,'transition');
end
Stx11verN=Stx11ver;
% Replace Parameters
for i=1:size(GemParam,1)-2
Stx11verN=strrep(Stx11verN,OldParam(i),VStruct(i));
end
% Create New Parameter File
NewParam=fopen('Stx11_New.dat','w');
fprintf(NewParam,strcat(Stx11verN));
fclose(NewParam);

% Check Solution Model File
if sum(ismember(VarSize{2}(1,:),"Wij"))==0
     % Update Solution Model
     [~,b]=find(ismember(VarSize{1}(1,:),"Wij"));
     WCom=join([Wij(:,1),Variables(sum(double(VarSize{1,1}(2,1:b-1)))+1:sum(double(VarSize{1,1}(2,1:b))))]);
     SolutionModel=fileread('Solution_Model.dat');
     OldW='W(' + string(strip(extractBetween(SolutionModel,'W(','d3'))) + 'd3';
     SMNew=SolutionModel;
     for i=1:length(OldW)
     SMNew=strrep(SMNew,OldW(i),WCom(i));
     end
     NewParam=fopen('Solution_Model_New.dat','w');
     fprintf(NewParam,strcat(SMNew));
     fclose(NewParam);
else
    SolutionModel=fileread('Solution_Model.dat');
    NewParam=fopen('Solution_Model_New.dat','w');
    fprintf(NewParam,strcat(SolutionModel));
    fclose(NewParam);
end
    
end