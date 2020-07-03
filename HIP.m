clear;clc;close all

xPrecision = [0.1,0.1,2*pi/20,0.1];
[x,y] = meshgrid(0:xPrecision(1):2);
load('V.mat');
Min = [0,0,0,0];Precision =[0.1; 0.1; 2*pi/20; 0.1];

alpha = 0:xPrecision(4):2*pi;
v = 0:xPrecision(3):2;
dt = 0.1;
ph = 10;
Tf = 10;

% u
w = .7*[1;1;1;1;1;1;1;-2;-2;-4];
w = 1.5*ones(10,1);
a = [5;5;10;5;2;2;2;-12;-12;-12];

uPrecision = [0.1,1];
Uw = -3:uPrecision(1):3;
Ua = -15:uPrecision(2):10;


plot(x,y,'.')

%% Dynamic
X0 = [0,0,0,0];
X(1,:)  = dynamic (X0,a(1),w(1),dt);
Xshow(1,:) = discrete(X(1,:),xPrecision);
for k=1:Tf
X (k+1,:) = dynamic (X(k,:),a(k),w(k),dt);
Xshow(k+1,:) = discrete(X(k+1,:),xPrecision);
end
hold on;scatter(Xshow(:,1),Xshow(:,2),'ro')
figure;

%% Goal (from PF)
% [x,y,alpha,v]
% G1 = [2,0.8,0,0];
% G2 = [0,2,0,0];
% G3 = [1,0,0,0];
G1 = [1.3000    1.7000    1.6000         0];
G2 = [1,1,0,0];
G = [G1;G2];

%% start
L=length(Uw)*length(Ua);
Beta = [0.1; 100];
Gamma = [0.9 ;0.99];

PBeta = (1/size(Beta,1)/size(G,1))*ones(size(Beta,1)*size(G,1),1)';
PGamma = (1/size(Gamma,1))*ones(size(Gamma,1),1)';

% weight of particles for goal position
WG = (1/size(G,1))*ones(size(G,1),1,1);

P=zeros(L,Tf,size(Gamma,1),size(Beta,1),size(G,1));
% showing results
PB(1,:) = PBeta;
PG(1,:) = PGamma;
PW(1,:) = WG; 

for k=1:Tf
    % prediction of u before observation
    for m=1:length(Beta)
        beta = Beta(m);
        for n=1:length(Gamma)
            gamma = Gamma(n);
            for i=1:size(G,1)
               % [Q,U] = evaluateQ(X(k,:),Uw,Ua,dt,G1,gamma);
                [Q,U] = evaluateQ2(phi,Min,Precision,X,Uw,Ua,dt,G(i,:),gamma);
                % i > m > n :MSB
                %^P(:,k,4*(n-1)+2*(m-1)+i) = evaluateP(Q,beta); %P(:,k,n,m,i)
                P(:,k,n,m,i) = evaluateP(Q,beta);
            end
        end
    end
    
    % update parameters after observing u
    u = [w(k),a(k)];
    ia = find(U(:,1)==w(k));
    ib = find(U(:,2)==a(k));
    for t=1:length(ia)
        aa=find(ib==ia(t));
        if(~isempty(aa))
            break;
        end
    end
    ind = ib(aa);
    
    
    % update gamma
    % Find the Q for the observed u (for debuging)
    %[Q1,U1] = evaluateQ2(phi,Min,Precision,X(k,:),u(1),u(2),dt,G1,Gamma(1));
    %[Q2,U2] = evaluateQ2(phi,Min,Precision,X(k,:),u(1),u(2),dt,G1,Gamma(2));
    
    % expected
%     PP = [PBeta(1)*WG(1) PBeta(1)*WG(2) PBeta(2)*WG(1) PBeta(2)*WG(2)];
%     Pu1 = sum([P(ind,k,1) P(ind,k,2) P(ind,k,3) P(ind,k,4)].*PP);
%     Pu2 = sum([P(ind,k,5) P(ind,k,6) P(ind,k,7) P(ind,k,8)].*PP);
%     PGamma(1) = PGamma(1)*Pu1/(PGamma(1)*Pu1+PGamma(2)*Pu2);
%     PGamma(2) = PGamma(2)*Pu2/(PGamma(1)*Pu1+PGamma(2)*Pu2);
    
    PG1 = P(ind,k,1,1,1)+P(ind,k,1,2,1)+P(ind,k,1,1,2)+P(ind,k,1,2,2);
    PG2 = P(ind,k,2,1,1)+P(ind,k,2,2,1)+P(ind,k,2,1,2)+P(ind,k,2,2,2);
    % filter role
    PGamma(1) = PGamma(1)*PG1;
    PGamma(2) = PGamma(2)*PG2;
    % normalization
    SUM = sum(PGamma);
    PGamma = PGamma/SUM;
    % save purpose only
    PG(k+1,1) = PGamma(1);
    PG(k+1,2) = PGamma(2);
    
    % update beta
    % Find the Q for the observed u
%     [Q1,U1] = evaluateQ2(phi,Min,Precision, X(k,:),u(1),u(2),dt,G1,Gamma(1));
%     [Q2,U2] = evaluateQ2(phi,Min,Precision, X(k,:),u(1),u(2),dt,G1,Gamma(2));
%     PP = [PGamma(1)*WG(1) PGamma(1)*WG(2) PGamma(2)*WG(1) PGamma(2)*WG(2)];
%     Pb1g1 = sum([P(ind,k,1) P(ind,k,5)].*[PGamma(1) PGamma(2)]*WG(1));
%     Pb1g2 = sum([P(ind,k,2) P(ind,k,6)].*[PGamma(1) PGamma(2)]*WG(2));
%     Pb2g1 = sum([P(ind,k,3) P(ind,k,7)].*[PGamma(1) PGamma(2)]*WG(1));
%     Pb2g2 = sum([P(ind,k,4) P(ind,k,8)].*[PGamma(1) PGamma(2)]*WG(2));
%      sm = sum(PBeta.*[Pb1g1 P1g2 Pb2g1 Pb2g2]);


    Pb1g1 = P(ind,k,1,1,1)+P(ind,k,2,1,1);
    Pb1g2 = P(ind,k,1,1,2)+P(ind,k,2,1,2);
    Pb2g1 = P(ind,k,1,2,1)+P(ind,k,2,2,1);
    Pb2g2 = P(ind,k,1,2,2)+P(ind,k,2,2,2);
    % filtering role
    PBeta(1) = PBeta(1)*Pb1g1;
    PBeta(2) = PBeta(2)*Pb1g2;
    PBeta(3) = PBeta(3)*Pb2g1;
    PBeta(4) = PBeta(4)*Pb2g2;
    % normalization
    sm = sum(PBeta);
    PBeta = PBeta/sm;
  
    
    % save purpose only
    PB(k+1,1) = PBeta(1);
    PB(k+1,2) = PBeta(2);
    PB(k+1,3) = PBeta(3);
    PB(k+1,4) = PBeta(4);
    
    % update g
    WG(1) = PBeta(1)+PBeta(3);
    WG(2) = PBeta(2)+PBeta(4);
    % normalize 
    sm = sum(WG);
    WG(1) = WG(1)/sm;
    WG(2) = WG(2)/sm;
    
    % save purpose only
    PW(k+1,1) = WG(1);
    PW(k+1,2) = WG(2);
    % new goal sampling from goal distribution
    
    % expected position in 1 step ahead (dynamic)
    % not biased :: learning rules of movement rather than dynamic
    Xp = dynamic(X(k,:),U(:,2),U(:,1),dt);
    Pt = ones(size(PG,2),size(PB,2),size(PW,2));
    Pt(1,:,:) = Pt(1,:,:)*PG(k+1,1);
    Pt(2,:,:) = Pt(2,:,:)*PG(k+1,2);
    Pt(:,1,:) = Pt(:,1,:)*PB(k+1,1);
    Pt(:,2,:) = Pt(:,2,:)*PB(k+1,2);
    Pt(:,3,:) = Pt(:,3,:)*PB(k+1,3);
    Pt(:,4,:) = Pt(:,4,:)*PB(k+1,4);
    Pt(:,:,1) = Pt(:,:,1)*PW(k+1,1);
    Pt(:,:,2) = Pt(:,:,2)*PW(k+1,2);
%     PXp = P(:,k,:,:,:).*Pt;
%     PXp =sum(sum(sum(PXp,3),4),5);
%     % mean
%     XP_mean = sum(Xp.*PXp,'all');
end

[~,i]=max(P(:,:,1));
ff=U(i,:);
subplot(211);plot(ff(:,1));
hold on; plot(w,'r');
subplot(212);plot(ff(:,2));
hold on; plot(a,'r');
figure;
[~,i]=max(P(:,:,2));
ff=U(i,:);
subplot(211);plot(ff(:,1));
hold on; plot(w,'r');
subplot(212);plot(ff(:,2));
hold on; plot(a,'r');
figure;
subplot(221);plot(PB(:,1))
title('Probability of Beta')
ylabel('Beta = 1, G= g1')
subplot(222);plot(PB(:,2))
ylabel('Beta = 1, G= g2')
subplot(223);plot(PB(:,3))
ylabel('Beta = 10, G= g1')
subplot(224);plot(PB(:,4))
ylabel('Beta = 10, G= g2')

figure;
subplot(211);plot(PG(:,1))
ylabel('gamma = 0.9')
title('Probability of Gamma')
subplot(212);plot(PG(:,2))
ylabel('gamma = 0.99')

figure;
subplot(211);plot(PW(:,1))
ylabel('Goal 1 ')
title('Goal Probability')
subplot(212);plot(PW(:,2))
ylabel('Goal 2 ')

