function [x,fval,exitflag,output,Cov,CovSPD,yminEst] = fminsearchbnd(fun,x0,LB,UB,options,Hess,varargin)
% FMINSEARCHBND: FMINSEARCH, but with bound constraints by transformation 
% Modified to Calculate Hessian Matrix ( Nelder & Mead 1965). DK 2021 
% Modified to Calculate Nearest SPD of Hessian Matrix.        DK 2021
% fminsearchbnd                            SPD - Adapted from NearestSPD
% Copyright (c) 2006, John D'Errico        Copyright (c) 2013, John D'Errico
% All rights reserved.                     All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
% 
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in
%       the documentation and/or other materials provided with the distribution
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.

% usage: x=FMINSEARCHBND(fun,x0)
% usage: x=FMINSEARCHBND(fun,x0,LB)
% usage: x=FMINSEARCHBND(fun,x0,LB,UB)
% usage: x=FMINSEARCHBND(fun,x0,LB,UB,options)
% usage: x=FMINSEARCHBND(fun,x0,LB,UB,options,p1,p2,...)
% usage: [x,fval,exitflag,output]=FMINSEARCHBND(fun,x0,...)
% 
% arguments:
%  fun, x0, options - see the help for FMINSEARCH
%
%  LB - lower bound vector or array, must be the same size as x0
%
%       If no lower bounds exist for one of the variables, then
%       supply -inf for that variable.
%
%       If no lower bounds at all, then LB may be left empty.
%
%       Variables may be fixed in value by setting the corresponding
%       lower and upper bounds to exactly the same value.
%
%  UB - upper bound vector or array, must be the same size as x0
%
%       If no upper bounds exist for one of the variables, then
%       supply +inf for that variable.
%
%       If no upper bounds at all, then UB may be left empty.
%
%       Variables may be fixed in value by setting the corresponding
%       lower and upper bounds to exactly the same value.
%
% Notes:
%
%  If options is supplied, then TolX will apply to the transformed
%  variables. All other FMINSEARCH parameters should be unaffected.
%
%  Variables which are constrained by both a lower and an upper
%  bound will use a sin transformation. Those constrained by
%  only a lower or an upper bound will use a quadratic
%  transformation, and unconstrained variables will be left alone.
%
%  Variables may be fixed by setting their respective bounds equal.
%  In this case, the problem will be reduced in size for FMINSEARCH.
%
%  The bounds are inclusive inequalities, which admit the
%  boundary values themselves, but will not permit ANY function
%  evaluations outside the bounds. These constraints are strictly
%  followed.
%
%  If your problem has an EXCLUSIVE (strict) constraint which will
%  not admit evaluation at the bound itself, then you must provide
%  a slightly offset bound. An example of this is a function which
%  contains the log of one of its parameters. If you constrain the
%  variable to have a lower bound of zero, then FMINSEARCHBND may
%  try to evaluate the function exactly at zero.
%
%
% Example usage:
% rosen = @(x) (1-x(1)).^2 + 105*(x(2)-x(1).^2).^2;
%
% fminsearch(rosen,[3 3])     % unconstrained
% ans =
%    1.0000    1.0000
%
% fminsearchbnd(rosen,[3 3],[2 2],[])     % constrained
% ans =
%    2.0000    4.0000
%
% See test_main.m for other examples of use.
%
%
% See also: fminsearch, fminspleas
%
%
% Author: John D'Errico
% E-mail: woodchips@rochester.rr.com
% Release: 4
% Release date: 7/23/06

% size checks
xsize = size(x0);
x0 = x0(:);
n=length(x0);

if (nargin<3) || isempty(LB)
  LB = repmat(-inf,n,1);
else
  LB = LB(:);
end
if (nargin<4) || isempty(UB)
  UB = repmat(inf,n,1);
else
  UB = UB(:);
end

if (n~=length(LB)) || (n~=length(UB))
  error 'x0 is incompatible in size with either LB or UB.'
end

% set default options if necessary
if (nargin<5) || isempty(options)
  options = optimset('fminsearch');
end

% stuff into a struct to pass around
params.args = varargin;
params.LB = LB;
params.UB = UB;
params.fun = fun;
params.n = n;
% note that the number of parameters may actually vary if 
% a user has chosen to fix one or more parameters
params.xsize = xsize;
params.OutputFcn = [];

% 0 --> unconstrained variable
% 1 --> lower bound only
% 2 --> upper bound only
% 3 --> dual finite bounds
% 4 --> fixed variable
params.BoundClass = zeros(n,1);
for i=1:n
  k = isfinite(LB(i)) + 2*isfinite(UB(i));
  params.BoundClass(i) = k;
  if (k==3) && (LB(i)==UB(i))
    params.BoundClass(i) = 4;
  end
end

% transform starting values into their unconstrained
% surrogates. Check for infeasible starting guesses.
x0u = x0;
k=1;
for i = 1:n
  switch params.BoundClass(i)
    case 1
      % lower bound only
      if x0(i)<=LB(i)
        % infeasible starting value. Use bound.
        x0u(k) = 0;
      else
        x0u(k) = sqrt(x0(i) - LB(i));
      end
      
      % increment k
      k=k+1;
    case 2
      % upper bound only
      if x0(i)>=UB(i)
        % infeasible starting value. use bound.
        x0u(k) = 0;
      else
        x0u(k) = sqrt(UB(i) - x0(i));
      end
      
      % increment k
      k=k+1;
    case 3
      % lower and upper bounds
      if x0(i)<=LB(i)
        % infeasible starting value
        x0u(k) = -pi/2;
      elseif x0(i)>=UB(i)
        % infeasible starting value
        x0u(k) = pi/2;
      else
        x0u(k) = 2*(x0(i) - LB(i))/(UB(i)-LB(i)) - 1;
        % shift by 2*pi to avoid problems at zero in fminsearch
        % otherwise, the initial simplex is vanishingly small
        x0u(k) = 2*pi+asin(max(-1,min(1,x0u(k))));
      end
      
      % increment k
      k=k+1;
    case 0
      % unconstrained variable. x0u(i) is set.
      x0u(k) = x0(i);
      
      % increment k
      k=k+1;
    case 4
      % fixed variable. drop it before fminsearch sees it.
      % k is not incremented for this variable.
  end
  
end
% if any of the unknowns were fixed, then we need to shorten
% x0u now.
if k<=n
  x0u(k:n) = [];
end

% were all the variables fixed?
if isempty(x0u)
  % All variables were fixed. quit immediately, setting the
  % appropriate parameters, then return.
  
  % undo the variable transformations into the original space
  x = xtransform(x0u,params);
  
  % final reshape
  x = reshape(x,xsize);
  
  % stuff fval with the final value
  fval = feval(params.fun,x,params.args{:});
  
  % fminsearchbnd was not called
  exitflag = 0;
  
  output.iterations = 0;
  output.funcCount = 1;
  output.algorithm = 'fminsearch';
  output.message = 'All variables were held fixed by the applied bounds';
  
  % return with no call at all to fminsearch
  return
end

% Check for an outputfcn. If there is any, then substitute my
% own wrapper function.
if ~isempty(options.OutputFcn)
  params.OutputFcn = options.OutputFcn;
  options.OutputFcn = @outfun_wrapper;
end

% now we can call fminsearch, but with our own
% intra-objective function.

%% Addition - Extract Simplex
[xu,fval,exitflag,output,Simplex] = Fminsearch(@intrafun,x0u,options,params);
%% End

% undo the variable transformations into the original space
x = xtransform(xu,params);

% final reshape to make sure the result has the proper shape
x = reshape(x,xsize);

% Use a nested function as the OutputFcn wrapper
  function stop = outfun_wrapper(x,varargin)
    % we need to transform x first
    xtrans = xtransform(x,params);
    
    % then call the user supplied OutputFcn
    stop = params.OutputFcn(xtrans,varargin{1:(end-1)});
    
  end

%% Addition - Initiate Hessian Calculation 
% Check to see if Hessian should be calculated 
if Hess==1
% Extract Simplex properties     
fval=Simplex{1};
two2np1=Simplex{2};
vu=Simplex{3};
fv=Simplex{5};
varargin=Simplex{6};
funfcn=Simplex{7};
TolV=options.TolX;
% Transform Variables back to original values 
for i = 1 : size(vu,2)
    vt=xtransform(vu(:,i),params);
    v(:,i)=reshape(vt,xsize);
    v(:,i)=v(:,i);
end

% Calculate Hessian, Psuedo-Hessian and x at minimum 
[Cov,CovSPD,x,yminEst]=Quad(fval,two2np1,v,fv,fun,TolV,UB); 

else 
    Cov=[];
    CovSPD=[];
    yminEst=x;
end

%% End 
end % mainline end

% ======================================
% ========= begin subfunctions =========
% ======================================
function [fval,xtrans] = intrafun(x,params)
% transform variables, then call original function

% transform
xtrans = xtransform(x,params);
% and call fun
fval = feval(params.fun,reshape(xtrans,params.xsize),params.args{:});

end % sub function intrafun end

% ======================================
function xtrans = xtransform(x,params)
% converts unconstrained variables into their original domains

xtrans = zeros(params.xsize);
% k allows some variables to be fixed, thus dropped from the
% optimization.
k=1;
for i = 1:params.n
  switch params.BoundClass(i)
    case 1
      % lower bound only
      xtrans(i) = params.LB(i) + x(k).^2;
      
      k=k+1;
    case 2
      % upper bound only
      xtrans(i) = params.UB(i) - x(k).^2;
      
      k=k+1;
    case 3
      % lower and upper bounds
      xtrans(i) = (sin(x(k))+1)/2;
      xtrans(i) = xtrans(i)*(params.UB(i) - params.LB(i)) + params.LB(i);
      % just in case of any floating point problems
      xtrans(i) = max(params.LB(i),min(params.UB(i),xtrans(i)));
      
      k=k+1;
    case 4
      % fixed variable, bounds are equal, set it at either bound
      xtrans(i) = params.LB(i);
    case 0
      % unconstrained variable.
      xtrans(i) = x(k);
      
      k=k+1;
  end
end

end % sub function xtransform end

%% Addition Hessian Calculator 
function[Cov,CovSPD,VMin,fval] = Quad(fval,two2np1,v,fv,fun,TolV,UB)    
% Quadratic surface fit based on Nelder and Mead 1965 
% Ensure df > rounding errors
% If change in Simplex points is 0
if sum(any(v(:,1)-v(:,2:end)==0))>1
    warning('Simplex converged to v1. Recalculate final simplex at lower resolution')
    %Cov=[];
    %CovSPD=[];
    %yminEst=[];
    %fv=[];
    TolReset=TolV/1e3;
    v(:,2:end)=v(:,1) + (TolReset+TolReset).*rand(length(v(:,1)),length(v(:,1))) - TolReset;
%     return 
end

% Expand final Simplex, check for local minimum and readjust Simplex. 
vp=v;
simp=10^-16;
TOL=TolV*2;
TolReset=TolV/1e3;
for i = 2:two2np1(end)
while fv(i)-fv(1)<simp && fv(i)-fv(1)>=0 || sum(fv(1,i)==fv(1,1:i-1))>length(v)*.001
    vp(:,i)=vp(:,i)+vp(:,i)-vp(:,1);
%   vp(vp(:,1)-vp(:,i)>TOL,i)=vp(vp(:,1)-vp(:,i)>TOL,1)-TOL;
%   vp(vp(:,1)-vp(:,i)<-TOL,i)=vp(vp(:,1)-vp(:,i)<-TOL,1)+TOL;
    fv(i)=fun(vp(:,i));
    fv(isnan(fv))=fv(1);
    if sum((vp(:,1)-vp(:,i)<=-TOL)) + sum(vp(:,1)-vp(:,i)>=TOL) == length(vp(:,1))
        if any(sum(vp(:,i)==vp(:,1:i-1))==length(vp(:,1))) || sum(fv(1,i)==fv(1,1:i-1))>length(v)*.01
            vp(:,i)=vp(:,1) + (TolReset+TolReset).*rand(length(vp(:,1)),1) - TolReset;
        else
        break
        end
    end

    if any(abs(vp(:,i))>UB(1))
        vp(:,i)=vp(:,1) + (TolReset+TolReset).*rand(length(vp(:,1)),1) - TolReset;
    end
   
end
end

v =[vp(:,(fv==min(fv))),vp(:,(fv ~=min(fv)))];
fval=min(fv);
%%
% Calculate fv at additional NAP points 
k=0;
Pij=zeros(size(v,1),(((size(v,1)+1)*(size(v,1)+2))/2)-length(v));
aval=zeros(size(v,1),size(v,1));

for i = 1:two2np1(end)
    for j = i:two2np1(end)
        if i~=j
            k=k+1;
            Pij(:,k)=(v(:,i) + v(:,j)) .* 0.5;
            [aval(i,j-1)] =  fun(Pij(:,k));
        end
    end
end

% Kmax=size(Pij,2);
% while min(aval(aval>0)) < fval
% minval =  min(aval(aval>0));  
% [in,dex]=find((aval==minval));    
% PS = 0;
% for i = 1:size(v,2)-in
%     PS = PS+i;
% end
% PK=Kmax-(PS-dex+in-1);
% fv(:,1)=fun(Pij(:,PK));
% v(:,1)=Pij(:,PK);
% k=0;
% i=1;
% for j = 2:two2np1(end)
% k=k+1;
% Pij(:,k)=(v(:,i) + v(:,j)) .* 0.5;
% [aval(i,j-1)] =  fun(Pij(:,k)); 
% end
% fval=fv(:,1);
% end

aval(isnan(triu(aval)))=max(max(aval))*10;
    
% Fill triangular Matrix
y0ij=aval;
y0ij = (y0ij+y0ij') - eye(size(y0ij,1)).*diag(y0ij);
a0=min(fv);
yi=fv(2:end);

Pi=v(:,2:end);
P0=v(:,1);

% Solve for y = a0 + 2a'x + x'Bx,
for i=1:two2np1(end)-1
ai(i,1)= (2*y0ij(1,i))-((yi(i)+3*a0)/2);
end

for i =1:two2np1(end)-1
bii(i,i) = 2*(yi(i)+a0 - 2*y0ij(1,i));
end

for i = 1:two2np1(end)-1
    for j = i:two2np1(end)-1
if i~=j
bij(i,j) = 2*(y0ij(i+1,j) +a0 - y0ij(1,i) - y0ij(1,j));
end
    end
end

% Fill B Matrix 
B=[bij;zeros(1,length(bij))];
B=B+bii;
B = (B+B') - eye(size(B,1)).*diag(B);

Q=Pi-P0;

xminEst=-(pinv(B)*ai);
PminEst=P0-(Q*(pinv(B)*ai));
yminEst=a0+(ai'*(pinv(B)*ai));

Cov=0.5*Q*pinv(B)*Q';
CovSPD=SPD(Cov);
% Calculate Hessian 
% pinv gives pseudo-inverse if inverse of Matrix isn't possible 
Hessian=2*(pinv(Q)'*B*pinv(Q));
%ensure symmetry
Hessian=(Hessian+Hessian')/2;

% If Hessian is not SPD 
HessianSPD=SPD(Hessian);

% v at min  
VMin=v(:,1);
end

function [SPDMatrix] = SPD(H)
% SPD - the nearest symmetric positive definite Matrix to the Hessian 
% Highman 1988
% Calculate singular value decomposition 
[U,S,V] = svd(H);
M = V*S*V';
% Calculate nearest SPD
SPDMatrix = (H+M)/2;
% Make symmetric
SPDMatrix = (SPDMatrix + SPDMatrix')/2;
% test for positive definiteness 
Check = 1;
n = 0;
while Check ~= 0 && n<1e6
  [~,Check] = chol(SPDMatrix);
  n = n + 1;
% Correct for floating point error
  if Check ~= 0
    Eig = min(eig(SPDMatrix));
    if Eig==0
        Eig=1e-16;
    end
    SPDMatrix = SPDMatrix +  (eye(size(H))*(-Eig*n.^2 + eps(Eig)));
  end
end
end
%% End


