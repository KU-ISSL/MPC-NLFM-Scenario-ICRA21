function plot_control(x,plan,constraints,ref,sys,car1,video)
%PLOT_CONTROL Plot the control trajectory from the current state
%   save: (0,1)

nx = size(x,1);
H = 9;
track_width = 1;
xplan = zeros(nx,H+1);

figure(1)

%plot goal point
plot(ref(1),ref(2),'k*','MarkerSize',7,'LineWidth',0.9) 
hold on

%compute constraints
xplot_const1 = linspace(-1,71,141);
plot(xplot_const1,3*sin(0.2*xplot_const1)-track_width,'k')
plot(xplot_const1,3*sin(0.2*xplot_const1)+track_width,'k')
%include active_constraints
plot(constraints(1,:),constraints(2,:),'k','LineWidth',1.5);
plot(constraints(1,:),constraints(3,:),'k','LineWidth',1.5);


%compute planned trajectory
xplan(:,1) = x;
for j = 2:H+1
    xplan(:,j) = sys(xplan(:,j-1),plan(j-1,:)',zeros(nx,1));
end
plot(xplan(1,:),xplan(2,:),'r-')

%plot estimated current state
car1 = car1.update(x(1), x(2), x(3), x(4));
plot(car1);

hold off
axis([-1 71 -5 5])
if nargin > 6
    frame = getframe(gcf);
    writeVideo(video,frame);
end
pause(0.5)

end

