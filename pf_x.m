function [xestsir,stdsir,xpartires,xpartires_1step]=pf_x(meas,xpartiold,U,prob)
%%
% Nparti: number of particles
% meas: meas(1), meas(2) are the measured x,y position
% xpartiold: is the Nparti of previous particles
% u: control signal u(1)=w , u(2)=a

%%
global L;

% initilization

Nparti=L;
Nparti1=1/Nparti;
stdmeas=(0.01); % standard deviation of the measurement errors

stdmodel1=0.001; % standard deviation of the evolution model
stdmodel2=0.001; % standard deviation of the evolution model
stdmodel3=0.001; % standard deviation of the evolution model
stdmodel4=0.001; % standard deviation of the evolution model


dt=0.1;


% initial condition
xestsir=[0;0;0;0];stdsir=[1;1;1;1];
xpartinew=zeros(4,Nparti);
x_onestep=zeros(4,Nparti);
wparti=zeros(1,Nparti);
wpartin=zeros(1,Nparti);
xpartires=zeros(4,Nparti);
xpartires_1step=zeros(4,Nparti);
wpartires=zeros(1,Nparti);

%%%%%%%%%%%%%%%%%%%%%%%%
% Advancing Particles %
%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:Nparti
%%%%%%%%%%%%%%%%%%%%
% 1 - PREDICTION %
%%%%%%%%%%%%%%%%%%%%

F=[1 0 0 dt*cos(xpartiold(3,i));0 1 0 dt*sin(xpartiold(3,i)); 0 0 1 0; 0 0 0 1];

%F=[1 dt 0;0 1 0;0 0 1];
% two step prediction
x_onestep(:,i)=F*xpartiold(:,i)+[0;0;dt*U(i,1);dt*U(i,2)];
xpartinew(:,i) = F*x_onestep(:,i)+[0;0;dt*U(i,1);dt*U(i,2)]+randn*[stdmodel1;stdmodel2;stdmodel3;stdmodel4];

% else
%     xpartinew(:,i)=randn*[stdmodel;0];
%  end
  
%   fe(1,i)=Ke*sign(xpartinew(1,i))*(abs(xpartinew(1,i)))^(ne)+Be*sign(xpartinew(1,i))*(abs(xpartinew(1,i)))^(ne)*xpartinew(2,i)+randn*stdrw1;
%  if ( xpartinew(1,i)<0)
%      fe(1,i)=0;
%  end

      

% integrating the probability of each action with the stdmeas
stdmeas = stdmeas/prob(i);
 %%%%%%%%%%%%%%%%%%%%
% Weights %
%%%%%%%%%%%%%%%%%%%%
wparti(i)=exp(-((xpartiold(1,i)-meas(1))/stdmeas)^2-((xpartiold(2,i)-meas(2))/stdmeas)^2 ...
            -((xpartiold(3,i)-meas(3))/stdmeas)^2 -((xpartiold(4,i)-meas(4))/stdmeas)^2);


end


wtotal=sum(wparti);
wpartin=wparti./wtotal;
%%%%%%%%%%%%%%%%
% 2 - UPDATE %
%%%%%%%%%%%%%%%%
%xpartiw=xpartinew.*wparti;
%%%%%%%%%%%%%%%%%%%%%%
% Mean at time k+1 %
%%%%%%%%%%%%%%%%%%%%%%

%xestpf=xpartinew*wpartin';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Standard Deviation at time k+1 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3 - RESAMPLING %
%%%%%%%%%%%%%%%%%%%%%%
%1
cresa=zeros(Nparti,1);
uresa=zeros(Nparti,1);
cresa(1)=wpartin(1);
for i=2:Nparti
cresa(i)=cresa(i-1)+wpartin(i);
end
iresa=1;
uresa(1)=rand*Nparti1;
for j=1:Nparti
uresa(j)=uresa(1)+Nparti1*(j-1);
while uresa(j) > cresa(iresa)
iresa=iresa+1;
end
xpartires(:,j)=xpartinew(:,iresa);
xpartires_1step(:,j)=x_onestep(:,iresa);
wpartires(j)=Nparti1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Means and standard deviation %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xestsir(1)=mean(xpartires(1,:)); % x
xestsir(2)=mean(xpartires(2,:)); % y
xestsir(3)=mean(xpartires(3,:)); % theta
xestsir(4)=mean(xpartires(4,:)); % v



stdsir(1)=std(xpartires(1,:)); % x
stdsir(2)=std(xpartires(2,:)); % dx
stdsir(3)=std(xpartires(3,:)); % theta
stdsir(4)=std(xpartires(4,:)); % v


%fe_hat=xestsir(3);
%%%%%%%%%%%%%%%%%%%%%
% End of resampling %
%%%%%%%%%%%%%%%%%%%%%
%xpartiold=xpartinew;
end