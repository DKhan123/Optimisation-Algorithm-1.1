% ------------------------------------------------------------------
    function[Hessian] = Quad(fval,two2np1,v,sigma,fv,varargin,funfcn)    
     %     QUADRATIC SURFACE FITTING
     
     %     EXPAND THE FINAL SIMPLEX, IF NECESSARY, TO OVERCOME ROUNDING
     %     ERRORS.
     
simp = 1e-13; % Default tolerance test 

for i = two2np1
    df = abs(fv(i)-fval);
    while (df<simp)
       v(:,i) = v(:,1)+(1/sigma)*(v(:,i) - v(:,1));
       x(:,1) = v(:,i); fv(:,i) = funfcn(x,varargin{:});
       df = abs(fv(i)-fval);
    end
end

 %     FUNCTION VALUES ARE CALCULATED AT AN ADDITIONAL NAP POINTS.

k=0;
for i = 1:two2np1(end)
    for j = i:two2np1(end)
        if i~=j
            k=k+1;
            Pij(:,k)=(v(:,i) + v(:,j)) .* 0.5;
            [aval(i,j-1)] = funfcn(Pij(:,k),varargin{:});
        end
    end
end
% Fill Triangular Matrix
y0ij=aval;
y0ij = (y0ij+y0ij') - eye(size(y0ij,1)).*diag(y0ij);
a0=fv(1);
yi=fv(2:end);
Pi=v(:,2:end);
P0=v(:,1);

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

B=[bij;zeros(1,length(bij))];
B=B+bii;
B = (B+B') - eye(size(B,1)).*diag(B);

xminEst=-(B\ai);
Q=Pi-P0;
PminEst=P0-(Q*(B\ai));
yminEst=a0+(ai'*(B\ai));

%s2=yminEst/(length(two2np1)-1);

%2*MSE*Hessian if N<n, N-n==1

Hessian2=2*(inv(Q)'*B*inv(Q))*yminEst;
Hessian=(inv(Q)'*B*inv(Q));
    end