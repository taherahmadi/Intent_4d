% function [Qh,U ] = evaluateQ2(V,Min,Precision,X,Uw,Ua,dt,G,gamma)

close all;
clear;

dim = 4;
% load('v_whole.mat');
% Grid = xs_whole;
% phi = ttrValue_obs;

load('V.mat');
Grid = xs;
% ttrValue_obs = cat(3, phi, phi(:,:,1,:)); \\ do we need to do this? if so
% how to make grids for interpn?
% phi = ttrValue_obs;

dt = 0.2;

gamma = 1;
beta = 10;
horizon = 1;
uPrecision = [1, 1];
Uw = -1:uPrecision(1):1;
Ua = -1:uPrecision(2):1;

Tf =200;

G = [0,0,pi/2,0];

% theta = -pi+9*pi/15
for theta=pi/2:pi/15:pi
X = zeros(Tf,4);
Umax = zeros(Tf,2);

init_X = [5, -5, theta, 0];

X(1,:) = init_X;
for k=1:Tf
    [Q,U,Xtn] = evaluateQ2(phi,Grid,X(k,:),Uw,Ua,dt,G,gamma,horizon);
    [P] = evaluateP(Q,beta);
%     figure;
%     plot(P);
    [Pmax, Idx] = max(P,[],2);
    Umax(k,:) = U(Idx,:);
    X(k+1,:) = dynamic (X(k,:),Umax(k,2),Umax(k,1),dt);
    
    if((X(k+1,1)-G(1))^2 + (X(k+1,2)-G(2))^2 < 0.2)
%          &&(abs(X(k+1,3)-G(3)) < 2*pi/15)
%          &&(abs(X(k+1,4)-G(4)) < 0.5)
        break;
    end
    
    fig = figure(2);
    
    fig_axes = axes('Parent',fig);
    hold(fig_axes, 'on');
    hold on; scatter(X(k,1),X(k,2),'r','LineWidth',2);
    for i=1:size(Xtn,1)
        [delta_x, delta_y] = pol2cart(Xtn(i,3), Xtn(i,4));
        hold on; quiver(X(k,1),X(k,2),delta_x, delta_y,0,'linewidth',0.5, 'MaxHeadSize',0.2)
        grid on;
        set(fig_axes, 'XLim',[-8 8], 'YLim', [-8 8]);
    end
    [delta_x, delta_y] = pol2cart(Xtn(Idx,3), Xtn(Idx,4));
    hold on; scatter(X(k,1),X(k,2),'r','LineWidth',2);
    hold on; quiver(X(k,1),X(k,2),delta_x,delta_y,0,'linewidth',4, 'MaxHeadSize',1.5)
end

figure;
[delta_x_g, delta_y_g] = pol2cart(G(:,3),0.5);
quiver(G(:,1),G(:,2),delta_x_g,delta_y_g,0,'linewidth',1, 'MaxHeadSize',1);

velocity = X(k,4)
if (X(k,4))==0
    velocity = 0.5
end
[delta_x_g, delta_y_g] = pol2cart(X(k,3),velocity);
quiver(X(k,1),X(k,2),delta_x_g,delta_y_g,0,'linewidth',1, 'MaxHeadSize',1);

hold on; scatter(G(1),G(2),'b');
hold on; plot(X(:,1),X(:,2),'r','LineWidth',2);
title(sprintf('Theta:%0.3f, k=%d',theta,k));

end

G_synth = G;
U_synth = Umax;
X_synth = X;
save('synthetic_trajectory.mat', 'X_synth', 'U_synth', 'G_synth');

