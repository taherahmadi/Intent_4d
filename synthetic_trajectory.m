% function [Qh,U ] = evaluateQ2(V,Min,Precision,X,Uw,Ua,dt,G,gamma)

close all;
clear all;
load('V.mat');
V = phi;
dim = 4;
load('V.mat');
Grid = xs;

dt = 0.5;
G = [3, 3, 0, 0];
gamma = 1;
horizon = 1;
uPrecision = [0.25, 1];
Uw = -0.5:uPrecision(1):0.5;
Ua = [-1.5000   -0.5000    0   0.5000    1.5000];

Tf = 10;
X = zeros(Tf,4);
Umax = zeros(Tf,2);
init_X = [0;0;0;0];
X(1,:) = init_X;
for k=1:10
    [Q,U] = evaluateQ2(phi,Grid,X(k,:),Uw,Ua,dt,G,gamma,horizon);
    plot(Q)
     [Qmax, Idx] = max(Q,[],2);
     Umax(k,:) = U(Idx,:);
     X(k+1,:) = dynamic (X(k,:),Umax(k,1),Umax(k,2),dt);
end

