function[x,Hessian]=fmin(a)
x0=a;
  options = optimset('PlotFcns',@optimplotfval,'TolFun',1e-3,'TolX',1e-3,'MaxIter',100000);
  [x,fval,exitflag,output,Hessian] = fminsearchbnd(@Mini,x0,[],[],options,1);
  
   function z = Mini(x)
       z=(1-x(1)).^2 + 100*(x(2)-x(1).^2).^2;
   end
end

