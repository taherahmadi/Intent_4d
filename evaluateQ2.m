function [Qh,U ] = evaluateQ2(V,Min,Precision,X,Uw,Ua,dt,G,gamma)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
x = X(1);
y = X(2);
alpha = X(3);
v = X(4);
xo = x + v*cos(alpha)*dt;
yo = y + v*sin(alpha)*dt;
alphao = alpha+ Uw'*dt;
vo = v+ Ua'*dt ;

L=length(Uw)*length(Ua);
p = meshgrid(alphao,vo);
p=reshape(p',[L,1]);
q = meshgrid(vo,alphao);
q=reshape(q,[L,1]);

xoo = xo*ones(L,1) + q.*cos(p)*dt;
yoo = yo*ones(L,1) + q.*sin(p)*dt;
Xo = [xo*ones(L,1),yo*ones(L,1),p,q];
Xoo = [xoo,yoo,p,q];


u1 = meshgrid(Uw,Ua);
u1=reshape(u1',[L,1]);
u2 = meshgrid(Ua,Uw);
u2=reshape(u2,[L,1]);

U=[u1 u2];
%Qh = -0*vecnorm([u1,u2]')- 0*abs(q')- gamma*sqrt((Xoo(:,1)-G(1)).^2+(Xoo(:,2)-G(2)).^2)';

% Find the discrete value of dynamic
Xd = discrete(Xo,Precision');
Xd(Xd<0) = 0;
ind = (Xd - Min)./Precision'+1;
rdist = [1 1] - G(1:2); % considering (0,0) at origin
relative_dist = (rdist - Min(1:2))./Precision(1:2)';
Vs1 = V(int16(ind+[relative_dist 0 0]));
r = -dt;
Qh = r*ones(L,1) - gamma* Vs1; % converting TTR to RL (NTTR) 

end

