clear all; clc; close all;

% This ttr code is generated for DubinsCar3D with additive disturbances
% dx = v * cos(theta) + d
% dy = v * sin(theta) + d
% dtheta = omega
% dv = a
% Control: a and omega
% disturbance: d


global a_upper; global a_lower;
global omega_max;
global target_rad; global target_v;
global target_x; global target_y;
global target_phi;

a_upper = 10; a_lower = -15;
omega_max = 3;

% target setting
target_rad = 0.01;
target_x = 0;
target_y = 0;
target_v = 0;
target_phi = 0;

dim = 4;
Min = zeros(dim,1);
Max = zeros(dim,1);
Min(1) = -1;
Min(2) = -1;
Min(3) = 0;
Min(4) = 0;
Max(1) = 1;
Max(2) = 1;
Max(3) = 2*pi;
Max(4) = 2;

% dx = [0.1; 0.1; 2*pi/100];
% dimension will be 41x41x40
dx = [0.05; 0.05; 2*pi/20; 0.1];
% dx = [0.05; 0.05; 2*pi/40];
Max(3) = Max(3) - dx(3);

[xs N] = gridGeneration(dim, Min, Max, dx);

% initialization
phi = 10*ones(N(1),N(2),N(3), N(4));

% Target init value = 0
% phi(((xs(:,:,:,1) - target_x).^2 + (xs(:,:,:,2) - target_y).^2) <= target_rad^2) = 0;
% change to circle 
flag =     ((target_x - xs(:,:,:,:,1)).^2 + (target_y - xs(:,:,:,:,2)).^2) <= target_rad &...
       abs(target_phi - xs(:,:,:,:,3)) <= target_phi;
   %        abs(target_rad - xs(:,:,:,:,3)) <= 0 &...
phi(flag) = 0;

% Save coordinate
% xs_whole = cat(3, xs, xs(:,:,1,:)); (TODO)
% save("/Users/anjianli/Desktop/robotics/project/ttr-compute/view/data/coordinate.mat", "xs_whole")


% %LF sweeping
% mex mex_dubins_car_reach_avoid_4d.cpp; disp('dubins car 4d reach avoid mexing done!');
% % disp('mexing has already been done!');
% 
% numIter = 10000;
% TOL = 0.5;
% 
% startTime = cputime;
% tic;
% mex_dubins_car_reach_avoid_4d(phi,xs,dx,a_upper,a_lower,omega_max,numIter,TOL);
% toc;
% mex_dubins_car_reach_avoid_4d
% endTime = cputime;
% fprintf('Total execution time %g seconds\n', endTime - startTime);


%LF sweeping
mex mexLFsweep.cpp; disp('mexing done!');

% numIter = 10000;
numIter = 50;
%TOL = eps;
TOL = 0.1;

startTime = cputime;
tic;
% mexLFsweep(phi,xs,dx,alpha,beta,V1,V2,numIter,TOL);
mexLFsweep(phi,xs,dx,a_upper,a_lower,omega_max,numIter,TOL);
toc;


contour(xs(:,:,1,10,1),xs(:,:,1,10,2),phi(:,:,1,10), "ShowText", "on")
xlabel('X')
ylabel('Y')
save('V.mat','phi')
    