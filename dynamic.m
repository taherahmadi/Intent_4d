function [Xout] = dynamic(X,a,w,dt)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
x = X(1);
y = X(2);
alpha = X(3);
v = X(4);
xo = x + v*cos(alpha)*dt+0*a;
yo = y + v*sin(alpha)*dt+0*a;
alphao = alpha+ w*dt;
vo = v+ a*dt ;
Xout = [xo,yo,alphao,vo];
end

