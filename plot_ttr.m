close all;
D1 = Min(1):dx(1):Max(1)
D2 = Min(2):dx(2):Max(2)
D3 = Min(3):dx(3):Max(3)
D4 = Min(4):dx(4):Max(4)

%% X-Y
for theta=[1]
%     for v=[1,2,3,4,5,6,8,10,12,14,16,21,24]
      for v=1:25
        contour(xs(:,:,1,10,1),xs(:,:,1,10,2),phi(:,:,theta,v), "ShowText", "on")
        xlabel('X')
        ylabel('Y')
        title(sprintf('theta=%0.2f v=%0.2f',D3(theta),D4(v)))
        f = gcf;
        % Requires R2020a or later
        exportgraphics(f,sprintf('plot_ttr/x_y/theta_%d_v_%d.png',theta,v))
    end
end
%% X-V
% for y=[1,10,20,30,40,60,80]
%     for theta=[1,2,3,4,5,6,10,15,20]
%         contour(reshape(xs(:,1,1,:,1),81,21),reshape(xs(:,1,1,:,4),81,21),reshape(phi(:,y,theta,:),81,21), "ShowText", "on")
%         xlabel('X')
%         ylabel('V')
%         title(sprintf('y=%0.2f theta=%0.2f',D2(y),D3(theta)))
%         f = gcf;
%         exportgraphics(f,sprintf('plot_ttr/x_v/y_%d_theta_%d.png',y,theta))
%     end
% end
% %% Y-V
% for x=[1,10,20,30,40,60,80]
%     for theta=[1,2,3,4,5,6,10,15,20]
%         contour(reshape(xs(1,:,1,:,2),81,21),reshape(xs(1,:,1,:,4),81,21),reshape(phi(x,:,theta,:),81,21), "ShowText", "on")
%         xlabel('Y')
%         ylabel('V')
%         title(sprintf('x=%0.2f theta=%0.2f',D1(x), D3(theta)))
%         f = gcf;
%         exportgraphics(f,sprintf('plot_ttr/y_v/x_%d_theta_%d.png',x, theta))
%     end
% end
% %% Y-Theta
% for x=[1,10,20,30,40,60,80]
%     for v=[1,2,3,4,5,6,10,15,21]
%         contour(reshape(xs(1,:,:,10,2),81,20),reshape(xs(1,:,:,10,3),81,20),reshape(phi(x,:,:,v),81,20), "ShowText", "on")
%         xlabel('Y')
%         ylabel('Theta')
%         title(sprintf('x=%0.2f v=%0.2f',D1(x),D4(v)))
%         f = gcf;
%         exportgraphics(f,sprintf('plot_ttr/y_theta/x_%d_v_%d.png',x,v))
%     end
% end
% %% X-Theta
% for y=[1,10,20,30,40,60,80]
%     for v=[1,2,3,4,5,6,10,15,21]
%         contour(reshape(xs(:,1,:,1,1),81,20),reshape(xs(:,1,:,1,3),81,20),reshape(phi(:,y,:,v),81,20), "ShowText", "on")
%         xlabel('X')
%         ylabel('Theta')
%         title(sprintf('y=%0.2f v=%0.2f',D2(y),D4(v)))
%         f = gcf;
%         exportgraphics(f,sprintf('plot_ttr/x_theta/y_%d_v_%d.png',y,v))
%     end
% end
% 
% %% Theta-V
% 
% for x=[1,10,20,30,40,60,80]
%     for y=[1,10,20,30,40,60,80]
%         contour(reshape(xs(1,1,:,:,3),20,21),reshape(xs(1,1,:,:,4),20,21),reshape(phi(x,y,:,:),20,21), "ShowText", "on")
%         xlabel('Theta')
%         ylabel('V')
%         title(sprintf('x=%0.2f y=%0.2f',D1(x),D2(y)))
%         f = gcf;
%         exportgraphics(f,sprintf('plot_ttr/theta_v/x_%d_y_%d.png',x,y))
%     end
% end