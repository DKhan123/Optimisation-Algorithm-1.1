function [x fval history] = VariableOut(x0,OtherVar,Names,StablePhases,E,E_Error,LBP,UBP,Ti,Pi,Wij)
    history = [];
    options = optimset('PlotFcns',@optimplotfval,'TolFun',1e-4,'TolX',1e-5,'MaxIter',500000,'OutputFcn', @myoutput);
   [x fval] = fminsearchbnd(@Mini,x0,LBP,UBP,options);
   %[x fval] = fminsearch(@Mini,x0,options)
        
    function stop = myoutput(x,optimvalues,state);
        stop = false;
        if isequal(state,'iter')
          history = [history, x];
        end
    end
    
    function z = Mini(x)
    z =  Minimise(x,OtherVar,Names,StablePhases,E,E_Error,Ti,Pi,Wij);
    %z=(x(1)^2-x(2)+1)^2;
    end
end