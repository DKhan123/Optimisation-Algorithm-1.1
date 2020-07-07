% This script is used to calculate new mean parameter values using Bayesian Inference. 
% Inputs are the parameter array, min and max value for the free parameter and the free parameter of intrest. 

function[NewMean,NewStDev]=Bayes(MeanStix,MeanParam,StDevParam,Stdev)
    % Construct Prior
    Mean    = MeanStix;
    S_D     = Stdev;
    dx      = (3*S_D/100000);
    x       = [(Mean-3*S_D):dx:(Mean+3*S_D)];
    Prior   = normpdf(x,Mean,S_D);
    
    % Construct Likelihood Gaussian
    Likelihood  = normpdf(x,MeanParam,StDevParam);
    Likelihood  = trapz(Prior)*dx*Likelihood/(trapz(Likelihood)*dx);
    % plot(x,Likelihood)
    % hold on 
    % plot(x,Prior)
    % Condition that x2=x1 
    p_condition=trapz(Prior.*Likelihood)*dx; 
    Pr_Li_condition=Likelihood.*Prior/p_condition;
    NewMean=trapz(x.*Pr_Li_condition)*dx;
    NewStDev=sqrt(trapz((x-NewMean).^2.*Pr_Li_condition)*dx);
end
