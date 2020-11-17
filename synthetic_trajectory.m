% function [Qh,U ] = evaluateQ2(V,Min,Precision,X,Uw,Ua,dt,G,gamma)

close all;
clear;

dim = 4;
load('V.mat');
Grid = xs;
% ttrValue_obs = cat(3, phi, phi(:,:,1,:)); \\ do we need to do this? if so
% how to make grids for interpn?
% phi = ttrValue_obs;

dt = 0.5;
G = [-5,-5,-pi/2,0];
gamma = 1;
beta = 10;
horizon = 1;
uPrecision = [0.5, 0.5];
Uw = -1:uPrecision(1):1;
Ua = -1.5:uPrecision(2):1.5;

Tf =15;
X = zeros(Tf,4);
Umax = zeros(Tf,2);
init_X = [0, 0, 0, 0];
X(1,:) = init_X;

for k=1:Tf
    [Q,U] = evaluateQ2(phi,Grid,X(k,:),Uw,Ua,dt,G,gamma,horizon);
    [P] = evaluateP(Q,beta);
    figure;
    plot(P);
    [Pmax, Idx] = max(Q,[],2);
    Umax(k,:) = U(Idx,:);
    X(k+1,:) = dynamic (X(k,:),Umax(k,1),Umax(k,2),dt);
end

figure;
scatter(G(1),G(2),'b')
hold on; plot(X(:,1),X(:,2),'r','LineWidth',2);
G_synth = G;
U_synth = Umax;
X_synth = X;
save('synthetic_trajectory.mat', 'X_synth', 'U_synth', 'G_synth');

