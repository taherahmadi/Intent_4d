clear all; close all;

%% xy data
load('./data/vicon_hat_3_hat_3_translation.csv')
xy_data = smoothdata(vicon_hat_3_hat_3_translation(600:1000,2:3),1,'gaussian');
txy = [vicon_hat_3_hat_3_translation(600:1000,1), xy_data(1:end,1:2)];

%% calculate the heading angle:
myDiff = diff(xy_data);
derived_theta = atan2(myDiff(:,2), myDiff(:,1));

txyheading = [xy_data(1:end-1,:), derived_theta];
%% velocity from xy and theta from velocity
t = txy(:,1);
len = length(t);
for i=2:len-1
	dx_1 = xy_data(i,1)-xy_data(i-1,1);
    dy_1 = xy_data(i,2)-xy_data(i-1,2);
    dt_1 = t(i,1)-t(i-1,1);
    dvx_1 = dx_1/dt_1;
    dvy_1 = dy_1/dt_1;

    dx_2 = xy_data(i+1,1)-xy_data(i,1);
    dy_2 = xy_data(i+1,2)-xy_data(i,2);
    dt_2 = t(i+1,1)-t(i,1);
    dvx_2 = dx_2/dt_2;
    dvy_2 = dy_2/dt_2;
    
    dvx(i-1) = (dvx_1+dvx_2)/2;
    dvy(i-1) = (dvy_1+dvy_2)/2;
    
    xy_v(i-1,:) = xy_data(i,1:2);
    tv(i-1) = t(i,1); 
    v(i-1) = sqrt(dvx(i-1)^2+dvy(i-1)^2);
    thetav(i-1) = atan2(dvy(i-1),dvx(i-1));
end
%% angular velocity from thetav
len = length(tv);
for i=2:len-1
	dt_1 = tv(i)-tv(i-1);
    dtheta_1 = thetav(i) - thetav(i-1);
    w_1 = dtheta_1/dt_1;
    
    dt_2 = tv(i+1)-tv(i);
    dtheta_2 = thetav(i+1) - thetav(i);
    w_2 = dtheta_2/dt_2;
    
    xy_w(i-1,:) = xy_data(i+1,1:2);
    tw(i-1) = tv(i);
    w(i-1) = (w_1+w_2)/2;
end
%% acceleration from v
len = length(tv);
for i=2:len-1
    dt_1 = tv(i)-tv(i-1);
    dv_1 = v(i) - v(i-1);
    a_1 = dv_1/dt_1;
    
    dt_2 = tv(i+1)-tv(i);
    dv_2 = v(i+1) - v(i);
    a_2 = dv_2/dt_2;
    
    xy_a(i-1,:) = xy_data(i+1,1:2);
    ta(i-1) = tv(i);
    a(i-1) = (a_1+a_2)/2;
end
%% plot
[x_prime,y_prime] = pol2cart(thetav(2:end-1), a);
xy_prime = [x_prime',y_prime'];

figure;
hold on;
scatter(xy_a(:,1),xy_a(:,2),'.');
quiver(xy_a(:,1),xy_a(:,2),xy_prime(:,1),xy_prime(:,2), 'MaxHeadSize',0.001);
xlabel('x');
ylabel('y');
%