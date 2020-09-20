clear;clc;close all

xPrecision = [0.1, 0.1, 2*pi/20,0.1];
[x,y] = meshgrid(0:xPrecision(1):2);
load('V.mat');
Min = [0,0,0,0];Precision =[0.1; 0.1; 2*pi/20; 0.1];

alpha = 0:xPrecision(4):2*pi;
v = 0:xPrecision(3):2;
dt = 0.1;
ph = 10;
Tf = 10;

% u
% w = .7*[1;1;1;1;1;1;1;-2;-2;-4];
w = 1.5*ones(10,1);
a = [5;5;10;5;2;2;2;-12;-12;-12];


% reduced size of possible input
% uPrecision = [0.5,5];
uPrecision = [1,10];

Uw = -2:uPrecision(1):2;
Ua = -15:uPrecision(2):10;

% plot(x,y,'.')

% Dynamic
X0 = [0,0,0,0];
X(1,:)  = dynamic (X0,a(1),w(1),dt);
Xshow(1,:) = discrete(X(1,:),xPrecision);
for k=1:Tf-1
X (k+1,:) = dynamic (X(k,:),a(k),w(k),dt);
Xshow(k+1,:) = discrete(X(k+1,:),xPrecision);
end
% hold on;scatter(Xshow(:,1),Xshow(:,2),'ro')
% figure;

%% Goal (from PF)
% [x,y,alpha,v]
% G1 = [2,0.8,0,0];
% G2 = [0,2,0,0];
% G3 = [1,0,0,0];

G1 = [1.3, 1.7, 1.6, 0];
% G2 = [0.9, 0, 0, 0];
% G3 = [1, 0.9, 0, 0];
% G4 = [1.9, 1.9, 0, 0];
% G5 = [0.1, 0.1, 0, 0];
% G6 = [0.1, 1.9, 1.6, 0];
% G7 = [1.9, 0.1, 0, 0];
% G8 = [1, 1.9, 1.4, 1];
% G9 = [1.9, 0.3, 0, 1];

% G = [G1;G2;G3;G4;G5;G6;G7;G8;G9];
% G = [G1;G2];

n = 5;
xmin=0;
xmax=2;
G=xmin+rand(n,4)*(xmax-xmin);
G(1,:) = G1;
%% start
L=length(Uw)*length(Ua);
Beta = [0.1; 100];
Gamma = [0.009 ;0.99];

PBeta = (1/size(Beta,1))*ones(size(Beta,1),1)';
PGamma = (1/size(Gamma,1))*ones(size(Gamma,1),1)';

% weight of particles for goal position
PGoal = (1/size(G,1))*ones(size(G,1),1,1);

Pu=zeros(L,Tf,size(Gamma,1),size(Beta,1),size(G,1));
PX = zeros(Tf,L);
Xp = zeros(L,4,Tf);
% showing results
PBt(1,:) = PBeta;
PGt(1,:) = PGamma;
PWt(1,:) = PGoal; 

for k=1:Tf
    
    % prediction of u before observation
    for m=1:length(Beta)
        beta = Beta(m);
        for n=1:length(Gamma)
            gamma = Gamma(n);
            for i=1:size(G,1)
               % [Q,U] = evaluateQ(X(k,:),Uw,Ua,dt,G1,gamma);
                [Q,U] = evaluateQ2(phi,Min,Precision,X(k,:),Uw,Ua,dt,G(i,:),gamma);
                % i > m > n :MSB
                %^P(:,k,4*(n-1)+2*(m-1)+i) = evaluateP(Q,beta); %P(:,k,n,m,i)
                Pu(:,k,n,m,i) = evaluateP(Q,beta);
            end
        end
    end
    
    % update parameters after observing u
    u = [w(k),a(k)];
%     ia = find(U(:,1)==w(k));
%     ib = find(U(:,2)==a(k));
    % if the actual u is not in the control samples find nearest
    ia = find(min(abs(U(:,1)-w(k)))==abs(U(:,1)-w(k)));
    ib = find(min(abs(U(:,2)-a(k)))==abs(U(:,2)-a(k)));
    
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
 
    PG = zeros(length(Gamma),1);
    for gamma=1:length(Gamma)
        PG(gamma) = sum(sum(Pu(ind,k,gamma,:,:))); %+Pu(ind,k,gamma,2,1)+Pu(ind,k,gamma,1,2)+Pu(ind,k,gamma,2,2);
        % filter role
        PGamma(gamma) = PGamma(gamma)*PG(gamma);
    end
    
    % normalization
    SUM = sum(PGamma);
    PGamma = PGamma/SUM;
    
    for gamma=1:length(Gamma)
        % multi purpose only
        PGt(k+1,gamma) = PGamma(gamma);
    end
 
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
    
    Pbg = zeros(length(Beta) , size(G,1));    
    for goal=1:size(G,1)
        for beta=1:length(Beta)
            for gamma=1:length(Gamma)
                Pbg(beta, goal) = Pbg(beta, goal) + Pu(ind,k,gamma,beta,goal)*PGamma(gamma);
            end
        end
    end
%     Pb1g1 = P(ind,k,1,1,1)+P(ind,k,2,1,1);
%     Pb1g2 = P(ind,k,1,1,2)+P(ind,k,2,1,2);
%     Pb2g1 = P(ind,k,1,2,1)+P(ind,k,2,2,1);
%     Pb2g2 = P(ind,k,1,2,2)+P(ind,k,2,2,2);
    
    for beta=1:length(Beta)
        % filtering role
        PBeta(beta) = PBeta(beta)*sum(Pbg(beta,:));
    end
    
%     PBeta(1) = PBeta(1)*Pb1g1;
%     PBeta(2) = PBeta(2)*Pb1g2;
%     PBeta(3) = PBeta(3)*Pb2g1;
%     PBeta(4) = PBeta(4)*Pb2g2;
    
    % normalization
    sm = sum(PBeta);
    PBeta = PBeta/sm;
  
    
    % save purpose only
    for beta=1:length(Beta)
        PBt(k+1,beta) = PBeta(beta);
    end
    
    % update g
    for goal=1:size(PGoal,1)
        PGoal(goal) = PGoal(goal)*sum(Pbg(:, goal));
    end
    
    % normalize 
    sm = sum(PGoal);
    PGoal = PGoal./sm;
    
%     G(:,1:2) = G(:,1:2) .* PGoal 
%     G(:,1:2) = xmin+G(:,1:2)*(xmax-xmin)
    
%     v1 = std(G(:,1),PGoal)
%     m1 = sum(PGoal.*G(:,1))/sum(PGoal)
%     pd1 = makedist('normal','mu',m1,'sigma',v1);
%     
%     v2 = std(G(:,1),PGoal*100)
%     m2 = sum((PGoal*100).*G(:,2))/sum(PGoal*100)    
%     pd2 = makedist('normal','mu',m2,'sigma',v2);
%     
%     x_values = 1:9;
%     G(:,1) = pdf(pd1,x_values)
%     G(:,2) = pdf(pd2,x_values)

    
    % save purpose only
    for goal=1:size(PGoal,1)
        PWt(k+1,goal) = PGoal(goal);
    end
    % new goal sampling from goal distribution
    
    % expected position in 1 step ahead (dynamic)
    % not biased :: learning rules of movement rather than dynamic
    Xp(:,:,k) = twostepdynamic(X(k,:),U(:,2),U(:,1),dt);
    Pt = ones(size(Gamma,1),size(Beta,1),size(G,1));
    for gamma=1:size(Gamma,1)
        Pt(gamma,:,:) = Pt(gamma,:,:)*PGt(k+1,gamma);
    end    
    for beta=1:size(Beta,1)
        Pt(:,beta,:) = Pt(:,beta,:)*PBt(k+1,beta);
    end
    for goal=1:size(G,1)
        Pt(:,:,goal) = Pt(:,:,goal)*PWt(k+1,goal);
    end
    
    PXp = ones(size(Pu));
    for jj =1:size(Pu,1)
        PXp(jj,k,:,:,:) = Pt; 
    end
    PXp(:,k,:,:,:) = Pu(:,k,:,:,:).*PXp(:,k,:,:,:);
    PX(k,:) =sum(sum(sum(PXp(:,k,:,:,:),3),4),5);
%     % mean
%     XP_mean = sum(Xp.*PXp,'all');

    
end


[~,i]=max(Pu(:,:,1,1,1));
ff=U(i,:);
subplot(211);plot(ff(:,1));
hold on; plot(w,'r');
subplot(212);plot(ff(:,2));
hold on; plot(a,'r');
figure;

[~,i]=max(Pu(:,:,2));
ff=U(i,:);
subplot(211);plot(ff(:,1));
hold on; plot(w,'r');
subplot(212);plot(ff(:,2));
hold on; plot(a,'r');
figure;

subplot(221);plot(PBt(:,1))
title('Probability of Beta')
ylabel('Beta = 1')
subplot(222);plot(PBt(:,2))
ylabel('Beta = 10')

figure;
subplot(211);plot(PGt(:,1))
ylabel('gamma = 0.009')
title('Probability of Gamma')
subplot(212);plot(PGt(:,2))
ylabel('gamma = 0.99')

figure;
subplot(911);plot(PWt(:,1))
ylabel('Goal 1 ')
title('Goal Probability')
subplot(912);plot(PWt(:,2))
ylabel('Goal 2 ')
subplot(913);plot(PWt(:,3))
ylabel('Goal 3 ')
subplot(914);plot(PWt(:,4))
ylabel('Goal 4 ')
subplot(915);plot(PWt(:,5))
ylabel('Goal 5 ')

% subplot(916);plot(PWt(:,6))
% ylabel('Goal 6 ')
% subplot(917);plot(PWt(:,7))
% ylabel('Goal 7 ')
% subplot(918);plot(PWt(:,8))
% ylabel('Goal 8 ')
% subplot(919);plot(PWt(:,9))
% ylabel('Goal 9 ')


% h = figure;
% for t=1:Tf
%     hold on;
%     Xp1 = reshape(Xp(:,1,t),[length(Uw),length(Ua)]);
%     Xp2 = reshape(Xp(:,2,t),[length(Uw),length(Ua)]);
%     PXt = reshape(PX(t,:),[length(Uw),length(Ua)]);
%     surf(Xp1,Xp2,10*log10(PXt));
%   
%    newmap = jet;                    %starting map
% ncol = size(newmap,1);           %how big is it?
% zpos = 1 + floor(2/3 * ncol);    %2/3 of way through
% newmap(zpos,:) = [1 1 1];        %set that position to white
% colormap(newmap); 
% end
figure;
view(2);
hold on;plot(Xshow(:,1),Xshow(:,2),'r','LineWidth',2);
hold on;scatter(Xshow(:,1),Xshow(:,2),50,'ro','LineWidth',2);

hold on;scatter(G(1,1),G(1,2),50,'go','LineWidth',8);
text(G(1,1),G(1,2), 'G1', 'Fontsize', 10);

for i=2:length(G)
hold on;scatter(G(i,1),G(i,2),50,'yo','LineWidth',8);
text(G(i,1),G(i,2), 'G'+string(i), 'Fontsize', 10);    
end

grid;