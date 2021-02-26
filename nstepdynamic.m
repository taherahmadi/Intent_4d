function [Xout,Unstep,Xnstep] = nstepdynamic(horizon,X,a,w,dt,U, G,phi,Grid,gamma,beta)
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

xo2 = xo + vo.*cos(alphao)*dt+0*a;
yo2 = yo + vo.*sin(alphao)*dt+0*a;
alphao2 = alphao+ w*dt;
vo2 = vo+ a*dt ;
Xout = [xo2,yo2,alphao2,vo2];
%for n TODO!
Unstep = zeros(length(G),horizon);
Xnstep = zeros(length(G),horizon);

optimal = 1;
if horizon>1
    for i=1:length(G)
        for j=2:horizon
            [Q,U,Xtn] = evaluateQ2(phi,Grid,[xo,yo,alphao,vo],U,dt,G(i,:),gamma,1);
            [P] = evaluateP(Q,beta);

            if optimal == 1
            % select optimal action/determinisitc
                [Pmax, Idx] = max(P,[],2);
            else
            % sample from action probability distribution
                Idx = datasample(ndgrid(1:size(P,2)),1,'weights',P);
            end

            a = U(Idx,2);
            w = U(Idx,1);
            Unstep(i,j) = [a,w];
            
            xo2 = xo + vo.*cos(alphao)*dt+0*a;
            yo2 = yo + vo.*sin(alphao)*dt+0*a;
            alphao2 = alphao+ w*dt;
            vo2 = vo+ a*dt ;
%             Xout = [xo2,yo2,alphao2,vo2];
            Xnstep(i,j) =  [xo2,yo2,alphao2,vo2];
        end
    end
end


end

