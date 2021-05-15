%% ME451 Final Project
clear all; close all; clc; %#ok<CLALL>

%% Set Up

% Parameters
M = 0.5; % kg
m = 0.2; % kg
b = 0.1; % kg/s
g = 9.8; % m/s^2
l = 0.3; % m
I = (m*l^2)/3; % kgm^2 (uniform density)

% State-Space Equations

den = I*(M+m) + M*m*l^2;

A = [0,  1,                  0,              0;
     0,  -b*(I+m*l^2)/den,   (m*l)^2*g/den,  0;
     0,  0,                  0,              1;
     0,  -m*l*b/den,         m*g*l*(M+m)/den,0];
 
B = [0;(I+m*l^2)/den;0;m*l/den];
 
C = [1,0,0,0;
     0,0,1,0];
 
D = [0;0];

% Encoder is limiting factor w/ 150 Hz --> use 100 Hz to be safe
Ts = 0.01; % sample time in seconds

states = {'x' 'x_dot' 'phi' 'phi_dot'};
inputs = {'u'};
outputs = {'x'; 'phi'};

% Continuous-time system
sysc = ss(A,B,C,D,'statename',states,'inputname',inputs,'outputname',outputs);

% Discrete-time system (zero-order hold)
sysd = c2d(sysc,Ts,'zoh');

%% Continuous Time Conrollability and Observability

ctrb_matc = ctrb(sysc);
obsv_matc = obsv(sysc);

isControllablec = (size(ctrb_matc,2) == rank(ctrb_matc));
isObservablec = (size(obsv_matc,2) == rank(obsv_matc));

if isControllablec
    disp('This system is controllable!')
else
    disp('This system is not controllable. Give up :~(')
end

if isObservablec
    disp('This system is observable!')
else
    disp('This system is not observable. Give up :~(')
end

%% Continuous Time Pole Placement

Q = diag([100,10,1000,10]);
R = 100;

K = lqr(sysc.A,sysc.B,Q,R);

sysc_cl = ss(sysc.A-sysc.B*K,sysc.B,sysc.C,sysc.D,'statename',states,'inputname',inputs,'outputname',outputs);

%% Impulse Response of Continous Time System

% [y,t] = impulse(sysc_cl,0:Ts:10);
% 
% figure(1)
% plot(t,y,'Linewidth',1);
% legend_c = legend('$x(t)$','$\phi (t)$','location','best');
% set(legend_c,'Interpreter','latex');
% title('Pendulum Position and Angle','interpreter','latex')
% xlabel('Time (t)','interpreter','latex'); 
% ylabel('Magnitude','interpreter','latex'); 

%% Simulate Continuous Time System

u = @(x,ref) -K*(x-ref);
delta = @(x,ref) sysc.A*x + sysc.B*u(x,ref);

tspan1 = 0:0.01:5;
ref1 = [1;0;0;0];
x0 = [0;0;0;0];
[t1,x1] = ode45(@(t,x) delta(x,ref1),tspan1,x0);

tspan2 = 5.01:0.01:10;
ref2 = [2;0;0;0];
x0 = [1;0;0;0];
[t2,x2] = ode45(@(t,x) delta(x,ref2),tspan2,x0);

t = [t1;t2]; x = [x1;x2];
u1 = -K*(x1.'-ref1); u2 = -K*(x2.'-ref2); u = [u1, u2];

figure(2)
plot(t,x(:,1),t,x(:,3),'Linewidth',1);
legend_c = legend('$x(t)$','$\phi (t)$','location','best');
set(legend_c,'Interpreter','latex');
title('Pendulum Position and Angle','interpreter','latex')
xlabel('Time (t)','interpreter','latex'); 
ylabel('Magnitude','interpreter','latex'); 

figure(3)
plot(t,u,'Linewidth',1)
title('Control Response','interpreter','latex')
xlabel('Time ($t$)','interpreter','latex'); 
ylabel('Magnitude ($u(t)$)','interpreter','latex'); 

%% Create Video Simulation

figure(4)
filename = 'continuousTime.avi';
for index = 1:10:length(t)
    simulate_pendulum(x(index,:));
    video((index-1)/10 + 1) = getframe(gcf);
    hold off;
end

% Create video from plots
writerObj = VideoWriter(filename);
writerObj.FrameRate = 10;

open(writerObj);
for index=1:(length(t)-1)/10 + 1
    frame = video(index);    
    writeVideo(writerObj, frame);
end
close(writerObj);

%% Discrete Time Conrollability and Observability

ctrb_matd = ctrb(sysd);
obsv_matd = obsv(sysd);

isControllabled = (size(ctrb_matd,2) == rank(ctrb_matd));
isObservabled = (size(obsv_matd,2) == rank(obsv_matd));

if isControllabled
    disp('This system is controllable!')
else
    disp('This system is not controllable. Give up :~(')
end

if isObservabled
    disp('This system is observable!')
else
    disp('This system is not observable. Give up :~(')
end

%% Discrete LQR Pole Placement

Q = diag([100,10,1000,10]);
R = 100;

[K] = dlqr(sysd.A,sysd.B,Q,R);

sysd_cl = ss(sysd.A-sysd.B*K,sysd.B,sysd.C,sysd.D,Ts,'statename',states,'inputname',inputs,'outputname',outputs);

%% Simulate Discrete Time System

N = 1000;
ref = zeros(4,N); 
ref(1,1:N/2) = ones(1,N/2); ref(1,N/2+1:N) = 2*ones(1,N/2);

x0 = [0;0;0;0];
x = zeros(4,N); x(:,1) = x0;
y = zeros(2,N);
u = zeros(1,N);

for index = 1:N-1
    u(index) = -K*(x(:,index)-ref(:,index));
    x(:,index+1) = sysd.A*x(:,index) + sysd.B*u(index);
    y(:,index) = sysd.C*x(:,index);
end

figure(2)
plot(0:Ts:9.99,x(1,:),0:Ts:9.99,x(3,:),'Linewidth',1);
legend_c = legend('$x(t)$','$\phi (t)$','location','best');
set(legend_c,'Interpreter','latex');
title('Pendulum Position and Angle','interpreter','latex')
xlabel('Time (t)','interpreter','latex'); 
ylabel('Magnitude','interpreter','latex'); 

figure(3)
plot(0:Ts:9.99,u,'Linewidth',1)
title('Control Response','interpreter','latex')
xlabel('Time ($t$)','interpreter','latex'); 
ylabel('Magnitude ($u(t)$)','interpreter','latex'); 

%% Simulate Discrete Time System with Noise

rng(420);

N = 1000;
ref = zeros(4,N); 
ref(1,1:N/2) = ones(1,N/2); ref(1,N/2+1:N) = 2*ones(1,N/2);

Qv = 1e-3*eye(4); Qw = 1e-3*eye(2);
K_true = dare(sysd.A',sysd.C',Qv,Qw);

x_noisy = zeros(4,N); x_noisy(:,1) = x0;
y_noisy = zeros(2,N);
u_noisy = zeros(1,N);

v = randn(4,N); w = randn(2,N);

x_est = zeros(4,N); x_est(:,1) = x0;
x_pred = x0;
K_pred = zeros(4,4,N); K_pred(:,:,1) = eye(4);

for index = 1:N-1
    u_noisy(index) = -K*(x_noisy(:,index)-ref(:,index));
    x_noisy(:,index+1) = sysd.A*x_noisy(:,index) + sysd.B*u_noisy(index) + Qv*v(:,index);
    y_noisy(:,index) = sysd.C*x_noisy(:,index) + Qw*w(:,index);
    [x_est(:,index+1),x_pred,~,K_pred(:,:,index+1)] = ...
        covarianceKalman(x_pred,K_pred(:,:,index),y_noisy(:,index),u_noisy(:,index),...
        sysd_cl.A,sysd_cl.B,sysd_cl.C,Qv,Qw);
end

K_diff = zeros(N,1);
for index = 1:N
    K_diff(index) = norm(K_pred(:,:,index) - K_true);
end

%% Plot Discrete Time System

tspan = 0:Ts:10-Ts;

figure(5)

subplot(1,2,1)
plot(tspan,x(1,:),tspan,x_noisy(1,:),tspan,x_est(1,:),'LineWidth',1);
legend_kalman = legend('Ideal State','True State','Filtered State');
set(legend_kalman,'Interpreter','latex');
title('Position ($x(t)$)','interpreter','latex')
xlabel('Time','interpreter','latex'); 
ylabel('Magnitude','interpreter','latex');

subplot(1,2,2)
plot(tspan,x(3,:),tspan,x_noisy(3,:),tspan,x_est(3,:),'LineWidth',1);
legend_kalman = legend('Ideal State','True State','Filtered State');
set(legend_kalman,'Interpreter','latex');
title('Angle ($\phi(t)$)','interpreter','latex')
xlabel('Time','interpreter','latex'); 
ylabel('Magnitude','interpreter','latex');

figure(6)
subplot(1,2,1)
plot(1:N,vecnorm(x-x_est),1:N,vecnorm(x-x_noisy),'LineWidth',1)
legend_delta = legend('$\| x(t) - x_{noisy}(t)\|$','$\| x(t) - x_{est}(t)\|$');
set(legend_delta,'Interpreter','latex');
title('State Error Magnitude','interpreter','latex')
xlabel('Time','interpreter','latex'); 
ylabel('Magnitude','interpreter','latex');

subplot(1,2,2)
plot(1:N,K_diff,'LineWidth',1)
title('$\| K_{est}(t) - K_{true}\|$','interpreter','latex')
xlabel('Time','interpreter','latex'); 
ylabel('Magnitude','interpreter','latex');

%% Create Video Simulations

figure(7)
filename = 'discreteTime.avi';
for index = 1:10:length(tspan)
    simulate_pendulum(x(:,index));
    video((index-1)/10 + 1) = getframe(gcf);
    hold off;
end

% Create video from plots
writerObj = VideoWriter(filename);
writerObj.FrameRate = 10;

open(writerObj);
for index=1:(length(tspan)-1)/10 + 1
    frame = video(index);    
    writeVideo(writerObj, frame);
end
close(writerObj);

figure(8)
filename = 'noisydiscreteTime.avi';
for index = 1:10:length(tspan)
    simulate_pendulum(x_noisy(:,index));
    video((index-1)/10 + 1) = getframe(gcf);
    hold off;
end

% Create video from plots
writerObj = VideoWriter(filename);
writerObj.FrameRate = 10;

open(writerObj);
for index=1:(length(tspan)-1)/10 + 1
    frame = video(index);    
    writeVideo(writerObj, frame);
end
close(writerObj);

figure(9)
filename = 'filtereddiscreteTime.avi';
for index = 1:10:length(tspan)
    simulate_pendulum(x_est(:,index));
    video((index-1)/10 + 1) = getframe(gcf);
    hold off;
end

% Create video from plots
writerObj = VideoWriter(filename);
writerObj.FrameRate = 10;

open(writerObj);
for index=1:(length(tspan)-1)/10 + 1
    frame = video(index);    
    writeVideo(writerObj, frame);
end
close(writerObj);

%% Test Noise Limits

% Standardized setup
rng(420);

N = 500;
ref = zeros(4,1); 
x0 = [0;0;0;0];
v = randn(4,N); w = randn(2,N);

% Test position process noise
y_test = zeros(2,N);
ii = 0.42;
while sum(abs(y_test(2,:))>0.35)==0
    
    Qv = ii*diag([1,0,0,0]); Qw = zeros(2);
    x_test = x0;
    y_test = zeros(2,N);
    
    for index = 1:N-1
        u_test = -K*(x_test-ref);
        x_test = sysd.A*x_test + sysd.B*u_test + Qv*v(:,index);
        y_test(:,index) = sysd.C*x_test + Qw*w(:,index);
    end 
    
    disp(ii)
    ii = ii + 0.001;    
end

Qv_max_x = 0.427;

% Test linear velocity process noise
y_test = zeros(2,N);
ii = 0.17;
while sum(abs(y_test(2,:))>0.35)==0
    
    Qv = ii*diag([0,1,0,0]); Qw = zeros(2);
    x_test = x0;
    y_test = zeros(2,N);
    
    for index = 1:N-1
        u_test = -K*(x_test-ref);
        x_test = sysd.A*x_test + sysd.B*u_test + Qv*v(:,index);
        y_test(:,index) = sysd.C*x_test + Qw*w(:,index);
    end 
    
    disp(ii)
    ii = ii + 0.001;    
end

Qv_max_xdot = 0.178;

% Test angle process noise
y_test = zeros(2,N);
ii = 0.028;
while sum(abs(y_test(2,:))>0.35)==0
    
    Qv = ii*diag([0,0,1,0]); Qw = zeros(2);
    x_test = x0;
    y_test = zeros(2,N);
    
    for index = 1:N-1
        u_test = -K*(x_test-ref);
        x_test = sysd.A*x_test + sysd.B*u_test + Qv*v(:,index);
        y_test(:,index) = sysd.C*x_test + Qw*w(:,index);
    end 
    
    disp(ii)
    ii = ii + 0.0001;    
end

Qv_max_phi = 0.0288;

% Test angular velocity process noise
y_test = zeros(2,N);
ii = 0.34;
while sum(abs(y_test(2,:))>0.35)==0
    
    Qv = ii*diag([0,0,0,1]); Qw = zeros(2);
    x_test = x0;
    y_test = zeros(2,N);
    
    for index = 1:N-1
        u_test = -K*(x_test-ref);
        x_test = sysd.A*x_test + sysd.B*u_test + Qv*v(:,index);
        y_test(:,index) = sysd.C*x_test + Qw*w(:,index);
    end 
    
    disp(ii)
    ii = ii + 0.001;    
end

Qv_max_phidot = 0.342;

% Test angle observation noise
y_test = zeros(2,N);
ii = 0.11;
while sum(abs(y_test(2,:))>0.35)==0
    
    Qv = zeros(4); Qw = ii*diag([0,1]);
    x_test = x0;
    y_test = zeros(2,N);
    
    for index = 1:N-1
        u_test = -K*(x_test-ref);
        x_test = sysd.A*x_test + sysd.B*u_test + Qv*v(:,index);
        y_test(:,index) = sysd.C*x_test + Qw*w(:,index);
    end 
    
    disp(ii)
    ii = ii + 0.001;    
end

Qw_max_phi = 0.117;
