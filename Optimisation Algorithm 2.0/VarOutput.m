function [Cov, CovSPD, yminEst, x] = VarOutput(VarInputs,VarSize,ExperimentData,Wij,GemParam,BulkComp,H)
    x0=VarInputs(:,1);
    VarExponents=VarInputs(:,2);
    LBP=VarInputs(:,3);
    UBP=VarInputs(:,4);
    history = [];
    No_It=0;
    Var0=1e10;
    if H==0;
    options = optimset('TolFun',1e-3,'TolX',1e-3,'MaxIter',1e5,'OutputFcn', @myoutput,'PlotFcns',@optimplotfval);
    else
    options = optimset('TolFun',1e-3,'TolX',1e-1,'MaxIter',1e5,'OutputFcn', @myoutput,'PlotFcns',@optimplotfval);
    end
    if isempty(LBP)&&isempty(UBP)||abs(sum(LBP))==inf&&sum(UBP)==inf
%         x0=VarInputs(:,5);
%         VarInputs(:,6)=ones(length(x0),1);
%         VarInputs(:,5)=zeros(length(x0),1);
    end
    [x,fval,exitflag,output,Cov,CovSPD,yminEst] = fminsearchbnd(@Mini,x0,LBP,UBP,options,H);
    %[x,fval,exitflag,output,Hessian] = Fminsearch(@Mini,x0,options);

 
        
    function stop = myoutput(x,optimvalues,state);
        stop = false;
        if isequal(state,'iter')
          history = [history, x];
        end
    end
   
    function z = Mini(x)
    [z,Var_It] =  Minimum2(x,VarInputs,VarExponents,VarSize,ExperimentData,Wij,GemParam,BulkComp);
     
    if z<Var0
    No_It=round(No_It)+1;
    VarR(:,No_It)=[z;Var_It];
    hold on 
    scatter(No_It,VarR(1,end),'r')
    end

    end
end