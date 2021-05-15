clc; clear all; close all; %#ok<CLALL>

%% consts

%Assumed plant variables
M=0.5; %kg Mass of cart
m=0.2; %kg mass of pendulum
L=0.3; %m length of pendulum
g=9.81; %m/s^2 gravity
% w=.1; %m width of pendulum
I=1/3*L^2*m; %kg-m^2 mass moment of inertia
b=0.1; %   dampig coef


% robust Variation of variables
l_robus=.25:.05:.35; %m
m_robus=.18:.01:.22; %kg
M_robus=.4:.05:.55; %kg
I_robus=1/12*L^2*m_robus; %m^2-kg


t = (0:0.01:5); %seconds

%% TF's
s=tf('s');

%TF using variables for the assumed plant
% mapping pendulum angle to input force  TF
q=(M+m)*(I+m*L^2)-(m*L)^2;
tf_pen_num=m*L*s/q;
tf_pen_den=s^3 + b*(I+m*L^2)*s^2/q - (M+m)*m*g*L/q*s - b*m*g*L/q;
tf_pen=tf_pen_num/tf_pen_den;


%mapping cart position to input force  TF
tf_cart_num= ((I+m*L^2)*s^2 - m*g*L)/q;
tf_cart_den=s^4 + b*(I + m*L^2)*s^3/q - (M + m)*m*g*L*s^2/q - b*m*g*L*s/q;
tf_cart=tf_cart_num/tf_cart_den;

%TF1=cart, TF2= pendulum
sys_tf = [tf_cart ; tf_pen];


%angular frequency range
w=(10^-3:.001:100);


%loop over all permutations of variation and plot normilized diference
%between assumed and varied plant
for i=1:1:length(M_robus)
   for p=1:1:length(m_robus)
        for z=1:1:length(l_robus)
            % mapping pendulum angle to input force  TF
            q=(M_robus(i)+m_robus(p))*(I+m_robus(p)*l_robus(z)^2)-(m_robus(p)*l_robus(z))^2;
            tf_pen_num_robus=m_robus(p)*l_robus(z)*s/q;
            tf_pen_den_robus=s^3 + b*(I+m_robus(p)*l_robus(z)^2)*s^2/q - (M_robus(i)+m)*m_robus(p)*g*l_robus(z)/q*s - b*m_robus(p)*g*l_robus(z)/q;
            tf_pen_robus=tf_pen_num_robus/tf_pen_den_robus;
            
            %bode mag plot of normalized variation difference
            [mag,phase,wout] = bode((tf_pen_robus-tf_pen)/tf_pen,w);
            figure(1)
            semilogx(wout, squeeze(20*log10(mag)));
            hold on
        end
    end 
end

%% Adding weight function 

%set weight function parameters
%these parameters were modifie to get weight function slightly above all
%other bode mag plots
s_max=0.3162; %max value of weight functuion as w aproches inf
A_start=1.0000e-08; %starting value of weight functuion as w aproches inf
w0=.9; %bandwith
W=(s+w0*A_start)/(s/s_max+w0); %weight function

%add weight function to other figure 
[mag_w,phase_w,wout_w]=bode(W, w);
weight=semilogx(wout_w,squeeze(20*log10(mag_w)),'--','linewidth',3,'color','red');
legend(weight, 'Weight Function', 'Location', 'southeast')
title('Bode of the normilized differance between varied and expected plant')
xlabel('Frequency (rad/s)')
ylabel('Magnitude (dB)')
