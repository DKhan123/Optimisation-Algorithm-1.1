function [SPDMatrix] = SPD(H)
% SPD - the nearest symmetric positive definite Matrix to the Hessian 
% Highman 1988
% Calculate singular value decomposition 
[~,S,V] = svd(H);
M = V*S*V';
% Calculate nearest SPD
SPDMatrix = (H+M)/2;
% Make symmetric
SPDMatrix = (SPDMatrix + SPDMatrix')/2;
% test for positive definiteness 
Check = 1;
n = 0;
while Check ~= 0
  [~,Check] = chol(SPDMatrix);
  n = n + 1;
% Correct for floating point error
  if Check ~= 0
    Eig = min(eig(SPDMatrix));
    SPDMatrix = SPDMatrix +  (eye(size(H))*(-Eig*n.^2 + eps(Eig)));
  end
end
end