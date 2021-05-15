clc; clear all; close all; %#ok<CLALL>

%% consts

M=0.5; %kg Mass of cart
m=0.2; %kg mass of pendulum
L=0.3; %m length of pendulum
g=9.81; %m/s^2 gravity
% w=.1; %m width of pendulum
I=0.006;%1/12*L^2*m; %kg-m^2 mass moment of inertia
b=0.1; %   dampig coef
t = (0:0.01:5);

%% TF's
s=tf('s');

% mapping pendulum angle to input force  TF
q=(M+m)*(I+m*L^2)-(m*L)^2;
tf_pen_num=m*L*s/q;
tf_pen_den=s^3 + b*(I+m*L^2)*s^2/q - (M+m)*m*g*L/q*s - b*m*g*L/q;
tf_pen=tf_pen_num/tf_pen_den;


%mapping crt position to input force  TF
tf_cart_num= ((I+m*L^2)*s^2 - m*g*L)/q;
tf_cart_den=s^4 + b*(I + m*L^2)*s^3/q - (M + m)*m*g*L*s^2/q - b*m*g*L*s/q;
tf_cart=tf_cart_num/tf_cart_den;

%TF1=cart, TF2= pendulum
sys_tf = [tf_cart ; tf_pen];


%% H_inf controlers

%set weight function parameters
s_max=1; %max value of weight functuion as w aproches inf
A_start=.001; %starting value of weight functuion as w aproches inf
w0=40; %bandwith
W_s_inv=(s+w0*A_start)/(s/s_max+w0); %weight function inverse
W_s=(s/s_max+w0)/(s+w0*A_start); %weight function

%h_inf controler
[K,CL,GAM,INFO] = mixsyn(sys_tf(2),W_s,[],W_s_inv);

%step of close loop control
figure
step(CL, t);

%create input step function 
u=zeros(1,length(t));
for index = ceil(length(t)/2):length(t)
    u(1,index)=1;
end

%create noisy input step function (gause noise)
u_noise =awgn(u,20, 'measured');

%sim step with and without noise
[y,time,x] = lsim(CL,u,t);
[y_noise,time,x] = lsim(CL,u_noise,t);

figure
plot(t,y)
xlabel('time')
ylabel('system output')
legend('cart', 'pendulum');

figure
plot(t,y_noise)
xlabel('time')
ylabel('system output')
legend('cart', 'pendulum');

figure
plot(t,u, t, u_noise)
xlabel('time')
ylabel('input signal')
legend('no noise', 'noise');


%Make video of pendulum step function 
% figure
% filename = 'discreteTime.avi';
% for index = 1:10:length(time)
%     simulate_pendulum(x(index,:));
%     video((index-1)/10 + 1) = getframe(gcf);
%     hold off;
% end
% 
% Create video from plots
% writerObj = VideoWriter(filename);
% writerObj.FrameRate = 10;
% 
% open(writerObj);
% for index=1:(length(time)-1)/10 + 1
%     frame = video(index);    
%     writeVideo(writerObj, frame);
% end
% close(writerObj);


figure
bodemag(CL(2), W_s_inv)
legend('CL Pendulum Angle Mag Response', 'Weight Function') 

figure
margin(CL(1))
figure
margin(CL(2))




%% pid

%%FAILED IMMC HEURISTIC 
% eps=0.1;
% 
% tau=65;
% wn=sqrt(g/L);
% zeta=b/2/wn;
% k=.02;
% theta=0.3;
% tau_c=theta;
% kc=2*tau*zeta/(k*(tau_c+theta));
% tau_i=2*tau*zeta;
% tau_d=20/3*tau*zeta;
% cont_filt=1/(eps*tau_d*s+1);
% %k_pid=kc*(1+1/(tau_i*s)+tau_d*s);

%zeigler nichols tuning 
%gain where starts to oscilate
K_u=200;

%assosiated period of osscilation 
T_u=.21; 

K_p=.6*K_u;
K_i=1.2*K_u/T_u;
K_d=0.075*K_u*T_u;

k_pid=pid(K_p, 1, K_d);

% PID from Zeigler testing (umich website)
%k_pid=pid(100, 1, 20);


%make CL function
pid_cl=k_pid*sys_tf(2)/(1+sys_tf(2)*k_pid);%feedback(sys_tf(2),k_pid); %sys_tf(2)/(1+sys_tf(2)*k_pid);
pid_cl_cart=k_pid*sys_tf(1)/(1+sys_tf(1)*k_pid);%feedback(1, sys_tf(2)*k_pid)*sys_tf(1);

PID_cl=[pid_cl_cart, pid_cl];

figure
impulse(PID_cl,t)
%xlim([0 2.5])
%ylim([-.2 .2])
title('System impulse resp')
% xlable('time')
% ylabel('angle')

% figure
% impulse(pid_cl_cart,t)
% xlim([0 5])
% title('cart position impulse resp')
% % xlable('time')
% % ylabel('angle')


%sim with noise and without
[y_pid,time,x] = lsim(PID_cl(2),u,t.');
[y_pid_cart,time,x] = lsim(PID_cl(1),u,t.');
[y_pid_noise,time,x] = lsim(PID_cl(2),u_noise,t.');
[y_pid_noise_cart,time,x] = lsim(PID_cl(1),u_noise,t.');

figure
plot(t,y_pid_cart, t,y_pid_noise_cart)
xlabel('time')
ylabel('system output')
legend('cart', 'cart with noise');
title('Step of Cart with and without noise')

figure
plot(t,y_pid, t,y_pid_noise)
xlabel('time')
ylabel('system output')
legend('pendulum', 'pendulum with noise');
title('Step of pendulum with and without noise')


figure
margin(pid_cl)
figure
margin(pid_cl_cart)


%% lqr

%import LQR continious
LQR_SS_file = load('sysc_cl.mat');
A_ss= LQR_SS_file.sysc_cl.A;
B_ss=LQR_SS_file.sysc_cl.B;
C_ss= LQR_SS_file.sysc_cl.C;
D_ss=LQR_SS_file.sysc_cl.D;

[a, b]=ss2tf(A_ss, B_ss,C_ss,D_ss,1);
lqr_cl=tf(a(2,:),b);
lqr_cl_cart=tf(a(1,:),b);

LQR_tf=[lqr_cl_cart, lqr_cl];

%import LQR discrete
LQR_SSd_file = load('sysd_cl.mat');
A_ssd= LQR_SSd_file.sysd_cl.A;
B_ssd=LQR_SSd_file.sysd_cl.B;
C_ssd= LQR_SSd_file.sysd_cl.C;
D_ssd=LQR_SSd_file.sysd_cl.D;

[a_d, b_d]=ss2tf(A_ssd, B_ssd,C_ssd,D_ssd,1);
lqr_cl_d=tf(a_d(2,:),b_d);
lqr_cl_d_cart=tf(a_d(1,:),b_d);

LQR_d_tf=[lqr_cl_d_cart, lqr_cl_d];



%step sym 
figure
lsim(LQR_SS_file.sysc_cl,u,t);
figure
lsim(LQR_SS_file.sysc_cl,u_noise,t);

%plots for LQR continious
figure()
impulse(LQR_SS_file.sysc_cl,t)


figure
margin(LQR_tf(1))
figure
margin(LQR_tf(2))

%Discrete 

%step w/ and w/out noise
figure
lsim(LQR_SSd_file.sysd_cl,u,t);
figure
lsim(LQR_SSd_file.sysd_cl,u_noise,t);


figure()
impulse(LQR_SSd_file.sysd_cl,t)


figure
margin(LQR_d_tf(1))
figure
margin(LQR_d_tf(2))
%% Robust
%set weight function parameters
s_max=0.3162; %max value of weight functuion as w aproches inf
A_start=1.0000e-08; %starting value of weight functuion as w aproches inf
w0=.9; %bandwith
W=(s+w0*A_start)/(s/s_max+w0); 

%loop gain
L_hinf=K*sys_tf(2);
L_PID=k_pid*sys_tf(2);

%complementary sensitivity function
T_hinf=L_hinf/(1+L_hinf);
T_PID=L_PID/(1+L_PID);




figure
bodemag(1/W)
hold on 
bodemag(T_hinf)
hold on 
bodemag(T_PID)
hold on 
bodemag(lqr_cl)
hold on 
bodemag(lqr_cl_d)
legend('1/W', 'T for H{\infty}', 'T for PID', 'LQR' , 'LQR Discrete')
title('Bode plot for Robustness')