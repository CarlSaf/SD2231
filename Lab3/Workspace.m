clear variables
clc, close all

lcn         =   'southeast';
mp          =   0.16;   %kg
cp          =   0.4;    %Ns/m
kp          =   6.32;   %N/m
Time        =   3;     %Sec
s           =   tf('s');
dp=0;
dd=cp;

cp          =   (2*sqrt(kp*mp));    %Z=1

w           =   sqrt(kp/mp);
z           =   cp/(2*sqrt(kp*mp));

sinefreq    =   w/2;      %Rad/sec    

Tf          =   (cp*s+kp)/(mp*s^2+cp*s+kp);
Tzw         =   (2*z*w*s + w^2)/(s^2 + 2*z*w*s + w^2);
Tzw_cp0     =   (w^2)/(s^2 + w^2);
% bode (Tzw, Tzw_cp0, Tf);

hp=0; hd=1.3; hi=0;
mp=0.16; ms=0.16;
cp=0.8; cs=0.05;
kp=6.32; ks=0.0632;
T=0;


%% Simulink
% z=1;
% sim('Lab3', Time)
% Time_sine=Sim_sine.time;  sine_value=Sim_sine.signals.values(:,1:2);
% Time_step=Sim_sine.time; step_value=Sim_step.signals.values(:,1:2);
% 
% figure, title('Response of Transfer func 1')
% subplot(2,1,1), plot(Time_sine,sine_value*100, 'linewidth', 2)
% hold on, xlabel('Time [s]'), ylabel('Displacement [cm]')
% subplot(2,1,2), plot(Time_step,step_value*100, 'linewidth', 2)
% hold on, xlabel('Time [s]'), ylabel('Displacement [cm]')
% 
% z=1.5;
% sim('Lab3', Time)
% Time_sine=Sim_sine.time; sine_value=Sim_sine.signals.values(:,1);
% Time_step=Sim_sine.time; step_value=Sim_step.signals.values(:,1);
% 
% title('Response of Transfer func 1')
% subplot(2,1,1), plot(Time_sine,sine_value*100, 'linewidth', 2)
% hold on
% subplot(2,1,2), plot(Time_step,step_value*100, 'linewidth', 2)
% hold on
% 
% z=0.5;
% sim('Lab3', Time)
% Time_sine=Sim_sine.time; sine_value=Sim_sine.signals.values(:,1);
% Time_step=Sim_sine.time; step_value=Sim_step.signals.values(:,1);
% 
% subplot(2,1,1), plot(Time_sine,sine_value*100, 'linewidth', 2)
% legend('Zp \zeta = 1','Zw - Road','Zp \zeta = 1.5','Zp \zeta = 0.5','Location','southwest')
% hold on
% subplot(2,1,2), plot(Time_step,step_value*100, 'linewidth', 2)
% legend('Zp \zeta = 1','Zw - Road','Zp \zeta = 1.5','Zp \zeta = 0.5','Location','northwest')
% hold on
% 
% 
% %% PSA
% v=50;
% w=sinefreq/2;
% omega=1e-2:1e2;
% sw=4.028e-7/(2.88e-4 + 0.68*w^2 + w^4);
% 
% figure, title('PSA of Transfer func 1')
% xlabel('Disturbance frequency, logarithmic scale'), ylabel('Gain, logarithmic scale')
% 
% z=1;
% Tzw         =   (2*z*w*s + w^2)/(s^2 + 2*z*w*s + w^2);
% g = freqresp(Tzw,omega);
% sp=((abs(g(1,:))).^2).*sw*(1/v);
% semilogy(omega, sp, 'linewidth', 2)
% hold on
% 
% z=1.5;
% Tzw         =   (2*z*w*s + w^2)/(s^2 + 2*z*w*s + w^2);
% g = freqresp(Tzw,omega);
% sp=((abs(g(1,:))).^2).*sw*(1/v);
% semilogy(omega, sp, 'linewidth', 2)
% hold on
% 
% z=0.5;
% Tzw         =   (2*z*w*s + w^2)/(s^2 + 2*z*w*s + w^2);
% g = freqresp(Tzw,omega);
% sp=((abs(g(1,:))).^2).*sw*(1/v);
% semilogy(omega, sp, 'linewidth', 2)
% 
% title('PSA of Transfer func 1 at speed of 50m/s')
% xlabel('Disturbance frequency [rad/s]'), ylabel('Gain')
% legend('Zp \zeta = 1','Zp \zeta = 1.5','Zp \zeta = 0.5')


%% PD
cp          =   0.4;
z           =   cp/(2*sqrt(kp*mp));
z           =   cp/(2*sqrt(kp*mp));

dp=0;
dd=1.15;
Tf_active=kp/(mp*s^2+dd*s+dp+kp);
figure, step(Tf, Tf_active)

% figure
% dd=0.5
% for i=0:100
%     dd=dd+0.05
%     Tfc=kp/(mp*s^2+dd*s+dp+kp);
%     step(Tf, Tfc)
%     pause(0.5)
% end
sim('Lab3', Time)
Time_sine=Sim_sine.time; sine_value=Sim_sine.signals.values;
Time_step=Sim_sine.time; step_value=Sim_step.signals.values;

figure, title('Response of Transfer function')
subplot(2,1,1), plot(Time_sine,sine_value(:,1:3)*100, 'linewidth', 2)
hold on, xlabel('Time [s]'), ylabel('Displacement [cm]')
subplot(2,1,2), plot(Time_step,step_value(:,1:3)*100, 'linewidth', 2)
hold on, xlabel('Time [s]'), ylabel('Displacement [cm]')
legend('Zw - Road','Zp - Passive','Zp - Active','Location','northwest')

figure, bodemag(Tf, Tf_active)


%%
close all
hp=0; hd=1.3; hi=0;
TfPID=kp/(mp*s^2+hd*s+kp+hp+hi/s);

% figure
% for i=0:100
%     hd=hd+0.05
%     TfPID=kp*s/(mp*s^3+hd*s^2+kp*s+hp*s+hi);
%     step(Tf, TfPID)
%     pause(0.5)
% end

figure, step(Tf, TfPID)
sim('Lab3', Time)
Time_sine=Sim_sine.time; sine_value=Sim_sine.signals.values;
Time_step=Sim_sine.time; step_value=Sim_step.signals.values;

figure, title('Response of Transfer function')
subplot(2,1,1), plot(Time_sine,sine_value(:,1:4)*100, 'linewidth', 2)
hold on, xlabel('Time [s]'), ylabel('Displacement [cm]')
subplot(2,1,2), plot(Time_step,step_value(:,1:4)*100, 'linewidth', 2)
hold on, xlabel('Time [s]'), ylabel('Displacement [cm]')
legend('Zw - Road','Zp - Passive','Zp - Active','Zp - PID','Location','northwest')

figure, bodemag(Tf, TfPID)

%% Skyhook
close all
T=1.25;
Tf_sky=kp/(mp*s^2+T*s+kp);

% figure
% for i=0:100
%     T=T+0.05
%     Tf_sky=kp/(mp*s^2+T*s+kp);
%     step(Tf, Tf_sky)
%     pause(0.5)
% end

figure, step(Tf, TfPID)
sim('Lab3', Time)
Time_sine=Sim_sine.time; sine_value=Sim_sine.signals.values;
Time_step=Sim_sine.time; step_value=Sim_step.signals.values;

figure, title('Response of Transfer function')
subplot(2,1,1), plot(Time_sine,sine_value(:,1:5)*100, 'linewidth', 2)
hold on, xlabel('Time [s]'), ylabel('Displacement [cm]')
subplot(2,1,2), plot(Time_step,step_value(:,1:5)*100, 'linewidth', 2)
hold on, xlabel('Time [s]'), ylabel('Displacement [cm]')
legend('Zw - Road','Zp - Passive','Zp - Active','Zp - PID','Zp - Skyhook','Location','northwest')

figure, bodemag(Tf, TfPID)

%% Summary
close all
% figure, step(Tf, Tf_active, TfPID, Tf_sky)
% figure, bodemag(Tf, Tf_active, TfPID, Tf_sky)
% 
% sim('Lab3', Time)
% Time_sine=Sim_sine.time; sine_value=Sim_sine.signals.values;
% Time_step=Sim_sine.time; step_value=Sim_step.signals.values;
% 
% figure, title('Response of Transfer function')
% subplot(2,1,1), plot(Time_sine,sine_value(:,1:5)*100, 'linewidth', 2)
% hold on, xlabel('Time [s]'), ylabel('Displacement [cm]')
% subplot(2,1,2), plot(Time_step,step_value(:,1:5)*100, 'linewidth', 2)
% hold on, xlabel('Time [s]'), ylabel('Displacement [cm]')
% legend('Zw - Road','Zp - Passive','Zp - Active','Zp - PID','Zp - Skyhook','Location','northwest')
% 
% 
% v=50;
% w=sinefreq/2;
% sw=4.028e-7/(2.88e-4 + 0.68*w^2 + w^4);
% omega=0:25;
% figure
% g = freqresp(Tf,omega);
% sp=((abs(g(1,:))).^2).*sw*(1/v);
% semilogy(omega, sp, 'linewidth', 2)
% hold on
% 
% g = freqresp(Tf_active,omega);
% sp=((abs(g(1,:))).^2).*sw*(1/v);
% semilogy(omega, sp, 'linewidth', 2)
% hold on
% 
% g = freqresp(TfPID,omega);
% sp=((abs(g(1,:))).^2).*sw*(1/v);
% semilogy(omega, sp, 'linewidth', 2)
% hold on
% 
% g = freqresp(Tf_sky,omega);
% sp=((abs(g(1,:))).^2).*sw*(1/v);
% semilogy(omega, sp, 'linewidth', 2)
% 
% title('PSA of Transfer func 1 at speed of 50m/s')
% xlabel('Disturbance frequency [rad/s]')
% legend('Zp - Passive','Zp - Active','Zp - PID','Zp - Skyhook','Location','northwest')

%% 2DOF
mp=0.16; ms=0.16;
cp=0.8; cs=0.05;
kp=6.32; ks=0.0632;

nom=cp*cs*s^2+(ks*cp+kp*cs)*s+kp*ks;
den=mp*ms*s^4+(cs*mp+cp*ms+cs*ms)*s^3+(ks*mp+cp*ms+kp*ms+ms*ks)*s^2+(kp*cs+ks*cp)*s+kp*ks;
Tf_2dof=nom/den;

bodemag(Tf, Tf_2dof)

v=50;
w=sinefreq/2;
sw=4.028e-7/(2.88e-4 + 0.68*w^2 + w^4);
omega=0:25;
% figure
% g = freqresp(Tf,omega);
% sp=((abs(g(1,:))).^2).*sw*(1/v);
% semilogy(omega, sp, 'linewidth', 2)
% hold on
% 
% g = freqresp(TfPID,omega);
% sp=((abs(g(1,:))).^2).*sw*(1/v);
% semilogy(omega, sp, 'linewidth', 2)
% hold on
% 
% g = freqresp(Tf_2dof,omega);
% sp=((abs(g(1,:))).^2).*sw*(1/v);
% semilogy(omega, sp, 'linewidth', 2)
% title('PSA of Transfer func 1 at speed of 50m/s')
% xlabel('Disturbance frequency [rad/s]')
% legend('Passive - 1DOF','PID - 1DOF','Zp - 2DOF')%,'Location','northwest')
% 
% 
% sim('Lab3', Time)
% Time_sine=Sim_sine.time; sine_value=Sim_sine.signals.values;
% Time_step=Sim_sine.time; step_value=Sim_step.signals.values;
% 
% figure, title('Response of Transfer function')
% subplot(2,1,1), plot(Time_sine,sine_value(:,1:6)*100, 'linewidth', 2)
% hold on, xlabel('Time [s]'), ylabel('Displacement [cm]')
% subplot(2,1,2), plot(Time_step,step_value(:,1:6)*100, 'linewidth', 2)
% hold on, xlabel('Time [s]'), ylabel('Displacement [cm]')
% legend('Zw - Road','Zp - Passive','Zp - Active','Zp - PID','Zp - Skyhook','Zp - 2DOF')%,'Location','northwest')

%% Skyhook
Time=25;
T=0.1450;
sim('Skyhook', Time)
Time_sine=Sim_sine.time; sine_value=Sim_sine.signals.values;
Time_step=Sim_sine.time; step_value=Sim_step.signals.values;

A=[0 1 0 0;-ks/ms 0 ks/ms 0;0 0 0 1;ks/mp 0 -(ks+kp)/mp -cp/mp];
B=[0 0 0;1/ms 0 0;0 0 0;-1/mp kp/mp cp/mp];
C=[1 0 0 0];
D=[0 0 0];


% T=0;
% for i=0:100
%     T=T+0.005
%     sim('Skyhook', Time)
%     pause(0.5)
% end

figure, title('Response of Transfer function')
subplot(2,1,1), plot(Time_sine,sine_value*100, 'linewidth', 2)
hold on, xlabel('Time [s]'), ylabel('Displacement [cm]')
subplot(2,1,2), plot(Time_step,step_value*100, 'linewidth', 2)
hold on, xlabel('Time [s]'), ylabel('Displacement [cm]')
legend('Zw - Road','Zp - Passive','Zp - Active','Zp - PID','Zp - Skyhook','Zp - 2DOF')%,'Location','northwest')

TF=tf(ss(A,B,C,D));
Tf_2DOF_skyhook=TF(2);

nom=cp*cs*s^2+(ks*cp+kp*cs)*s+kp*ks;
den=mp*ms*s^4+(cs*mp+cp*ms+cs*ms)*s^3+(ks*mp+cp*ms+kp*ms+ms*ks)*s^2+(kp*cs+ks*cp)*s+kp*ks;
Tf_2dof=nom/den;

v=50;
w=sinefreq/2;
sw=4.028e-7/(2.88e-4 + 0.68*w^2 + w^4);
omega=0:25;

figure
g = freqresp(Tf_2dof,omega);
sp=((abs(g(1,:))).^2).*sw*(1/v);
semilogy(omega, sp, 'linewidth', 2)
hold on

g = freqresp(Tf_2DOF_skyhook,omega);
sp=((abs(g(1,:))).^2).*sw*(1/v);
semilogy(omega, sp, 'linewidth', 2)

title('PSA of Transfer func at speed of 50m/s')
xlabel('Disturbance frequency [rad/s]')
legend('Passive - 2DOF','Active - 2DOF')%,'Location','northwest')

figure, bodemag(Tf_2DOF_skyhook, Tf_2dof)
