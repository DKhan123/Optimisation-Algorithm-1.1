function [x fval history] = Wout(W,Ti,Pi,StablePhases,E,E_Error,LBW,UBW)
    history = [];
    options = optimset('PlotFcns',@optimplotfval,'TolFun',1e-5,'TolX',1e-5,'MaxIter',500000,'OutputFcn', @myoutput);
   [x fval] = fminsearchbnd(@f,W,LBW,UBW,options);
        
    function stop = myoutput(x,optimvalues,state);
        stop = false;
        if isequal(state,'iter')
          history = [history, x];
        end
    end
    
    function z = f(W)
    z =  W_Params(W,Ti,Pi,StablePhases,E,E_Error);
    end
end