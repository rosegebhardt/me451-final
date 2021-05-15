function simulate_pendulum(state)

hold off

% State Information
pos = state(1);
angle = -state(3) + pi;

% Dimensions
W = 1;
H = 0.5;
rh = 0.18;

% Cart Location
y = 0;
pendx = pos + 2*sin(angle);
pendy = y - 2*cos(angle);

% Draw cart
plot([-10 10],[0 0],'w','LineWidth',2), hold on
rectangle('Position',[pos-W/2,y-H/2,W,H],'Curvature',.1,'FaceColor','k','LineWidth',1.5);
% Draw pendulum
plot([pos pendx],[y pendy],'Color',[0.75 0.5 0],'LineWidth',6);
% Draw hinge
rectangle('Position',[pos-rh/2,y-rh/2,rh,rh],'Curvature',1,'FaceColor',[0.7,0.5,1],'LineWidth',0.5);

axis([-4 4 -1 3]); axis equal
set(gcf,'Position',[100 100 1000 400])
drawnow