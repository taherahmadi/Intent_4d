function [Qh,U ] = evaluateQ2(V,Grid,X,Uw,Ua,dt,G,gamma,horizon)


Length_U=length(Uw)*length(Ua);


u1 = meshgrid(Ua,Uw);
u1=reshape(u1,[Length_U,1]);
u2 = meshgrid(Uw,Ua);
u2=reshape(u2',[Length_U,1]);

U=[u1,u2];


    
Xt1 = dynamic(X,u2,u1,dt);




% n-step
for i=2:horizon
   % HARD
   % call evaluateP to greedy choose the current best action
   % for time t+horizen-1
end

if horizon ==1
    Xtn = Xt1;
end
% apply the translation and rotation due to goal change
Xtn(:,3) = Xtn(:,3)-G(3);
% check rotation
Xtn_rot = [cos(G(3)) sin(G(3));
      -sin(G(3)) cos(G(3))]*([Xtn(:,1)-G(1) Xtn(:,2)-G(2)])';
Xtn(:,1) = (Xtn_rot(1,:))';
Xtn(:,2) = (Xtn_rot(2,:))';


thet = Xtn(:,3) ;
thet(thet<0)=thet(thet<0)+2*pi;
thet(thet>2*pi)=thet(thet>2*pi)-2*pi;
Xtn(:,3) = thet;

% discount V before (geometric series)

x1 = Grid(:,:,:,:,1);
x2 = Grid(:,:,:,:,2);
x3 = Grid(:,:,:,:,3);
x4 = Grid(:,:,:,:,4);


TTR = interpn(x1, x2, x3, x4, V, Xtn(:,1), Xtn(:,2), Xtn(:,3), Xtn(:,4), 'nearest', 0);

%figure; plot(Vs1)



r = -dt;
V = r*ones(Length_U,1) - gamma*TTR;

if gamma < 1
    Qh = (((1-gamma^horizon)/(1-gamma)).*V)'; 
else
    Qh = V'; 
end
% Qh = -0*vecnorm([u1,u2]')- 0*abs(q')- gamma*sqrt((Xoo(:,1)-G(1)).^2+(Xoo(:,2)-G(2)).^2)';

end

