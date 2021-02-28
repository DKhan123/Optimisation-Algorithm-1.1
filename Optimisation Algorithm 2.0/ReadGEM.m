%% Read GEM Results 
function [GemData,TPVal]=ReadGEM(Result)
START=["Phase Compositions (molar  proportions):","Molar Properties and Density:","Seismic Properties:"];
END=["Phase speciation (molar proportions):","Seismic Properties:","Bulk Composition:"];

%BulkStart="Other Bulk Properties:";
%BulkEnd="Chemical Potentials (J/mol):";

SplitLine=cell(3,1);
for i =1:length(START)
List=extractBetween(Result,START(i),END(i));
    for k=1:length(List)
StringList=string(List(k));
LineList=splitlines(StringList);
LineList=LineList (~ strcmp (LineList (:,1), "" ), :);
LineList(1)='Phase'+LineList(1);
if contains(LineList(1),'Poisson ratio')
   LineList(1)= strrep(LineList(1),'Poisson ratio','Poisson_ratio');
end

for j=1:length(LineList)
Properties{i,k}(j,:)=(regexp(string(LineList(j)),'[\w(.|-)]+','match'));
end
Properties{i,k}(2:end,2:end)=str2double(Properties{i,k}(2:end,2:end));
    end
GemData=Properties;
end

Start=["Stable phases at:"];
End=["Phase Compositions"];
List=extractBetween(Result,Start,End);
for k=1:length(List)
StringList=string(List(k));
LineList=splitlines(StringList);
LineList=strip(LineList (~ strcmp (LineList (:,1), "" ), :));

for j=1:length(LineList)
Expression=(regexp(string(LineList(j)),'[\w(.|-)]+','match'));
TPVal(k,j)=double(Expression(end));
end
end

end
