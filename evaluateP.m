function [P] = evaluateP(Q,beta)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
eps = 0.000001;
P = exp(beta*Q)/(sum(exp(beta*Q))+eps);
end

