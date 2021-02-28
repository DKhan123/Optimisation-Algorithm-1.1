%% Read GEM Results 
function [GemData]=ReadGEM(Result)
START=["Molar Properties and Density:"];
END=["Seismic Properties:"];
%BulkStart="Other Bulk Properties:";
%BulkEnd="Chemical Potentials (J/mol):";

SplitLine=cell(1,1);
List=extractBetween(Result,START,END);
    for k=1:length(List)
StringList=string(List(k));
LineList=splitlines(StringList);
LineList=LineList (~ strcmp (LineList (:,1), "" ), :);
Array=extractBetween(strip(LineList(2:end-1)),1," ");

Array1{k}=Array;
%Array= arrayfun (@ (row) vertcat (Temp {row ,:}), (1: size (Temp, 1)), 'Uniform' , 0);
%Array1{k}= Array;
end
GemData=Array1;
end