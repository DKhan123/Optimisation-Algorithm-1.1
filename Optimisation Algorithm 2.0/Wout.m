function [x fval history VarR] = Wout(W,Ti,Pi,StablePhases,E,E_Error,LBW,UBW)
    history = [];
    No_It=0;
    Var0=1e10;
    options = optimset('PlotFcns',@optimplotfval,'TolFun',1e-5,'TolX',1e-5,'MaxIter',500000,'OutputFcn', @myoutput);
   [x fval] = fminsearchbnd(@f,W,LBW,UBW,options);
        
    function stop = myoutput(x,optimvalues,state);
        stop = false;
        if isequal(state,'iter')
          history = [history, x];
        end
    end
    
    function z = f(W)
    [z,Var_It] =  W_Params(W,Ti,Pi,StablePhases,E,E_Error);
    if z<Var0
    No_It=round(No_It)+1;
    VarR(:,No_It)=[z;Var_It];
    %Var0=VarR(1);
    end
    end
end