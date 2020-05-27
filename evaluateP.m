function [P] = evaluateP(Q,beta)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
P = exp(beta*Q)/sum(exp(beta*Q));
end

