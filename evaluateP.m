function [P] = evaluateP(Q,beta)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
% sm = sum(exp(beta*Q));
% 
% if sm <eps
%     P = exp(beta*Q*0.01)/sum(exp(beta*Q*0.01));
% else
%     P = exp(beta*Q)/sm;
% end

Q_max = max(Q);
a = exp(beta*(Q-Q_max));
b = sum(a);

P = a/b;
end