clc, close all, clear variables
Time=5;
testHZ=1
sinfreq=testHZ*2*pi;

m           =   22e3;
J           =   7e5;
c           =   4e4;
k           =   6e5;
L           =   6;


Ab = [0 1 0 0;
    -2*k/m -2*c/m 0 0;
    0 0 0 1;
    0 0 -2*L^2*k/J -2*L^2*c/J];

Bb = [0 0 0 0;
    k/m c/m k/m c/m;
    0 0 0 0;
    -k*L/J -c*L/J  -k*L/J -c*L/J];

Cb = [1 0 0 0;
    0 0 1 0];

Db = [0 0 0 0;
    0 0 0 0];

% TF=tf(ss(A,B,C,D));
% Tf_bounce=TF(1,:);
% Tf_pitch=TF(2,:);
% damp(Tf_bounce);
% damp(Tf_pitch);

% bodemag(Tf_pitch)

% sim('Bounce',Time)
% Time_sine=Sim_sine.time; sine_value=Sim_sine.signals.values;
% Time_step=Sim_sine.time; step_value=Sim_step.signals.values;
% 
% Sine_bounce     = sine_value(:,1:2);
% Sine_pitch      = sine_value(:,3);
% Step_bounce     = step_value(:,1:2);
% Step_pitch      = step_value(:,3);
% 
% figure
% subplot(2,1,1), plot(Time_sine,Sine_bounce*100, 'linewidth', 2)
% hold on, xlabel('Time [s]'), ylabel('Displacement [cm]')
% title('Response of Transfer function at: ' + string(testHZ) + ' Hz')
% subplot(2,1,2), plot(Time_step,Step_bounce*100, 'linewidth', 2)
% hold on, xlabel('Time [s]'), ylabel('Displacement [cm]')
% legend('Zw - Road','Z - Bounce')
% 
% 
% figure
% subplot(2,1,1), plot(Time_sine,rad2deg(Sine_pitch), 'linewidth', 2)
% hold on, xlabel('Time [s]'), ylabel('Displacement [\circ]')
% title('Response of Transfer function at: ' + string(testHZ) + ' Hz')
% subplot(2,1,2), plot(Time_step,rad2deg(Step_pitch), 'linewidth', 2)
% hold on, xlabel('Time [s]'), ylabel('Displacement [\circ]')
% legend('\chi - Pitch')

%%
% close all
m = 22000; %kg
J = 700000; %kgm^2
k = 600000;
L = 6;
Cx=4e6;
Cz=24e3;
Time=5;

A = [0 1 0 0;
    -2*k/m 0 0 0;
    0 0 0 1
    0 0 (-2*k*L^2)/J 0];

B = [0 0 0 0;
    k/m k/m -1/m -1/m;
    0 0 0 0;
    -k*L/J k*L/J L/J -L/J];

C = [0 1 0 0;
    0 0 0 1];

D = [0 0 0 0;
    0 0 0 0];

% Cz=1e5
% for i=0:100
%     Cz=Cz+5000
%     sim('Active_skyhook',Time)
%     pause(0.2)
% end

% Cx=2e6
% for i=0:100
%     Cx=Cx+50000
%     sim('Active_skyhook',Time)
%     pause(0.2)
% end

% sim('Active_skyhook',Time)
% Time_sine=Sim_sine.time; sine_value=Sim_sine.signals.values;
% Time_step=Sim_sine.time; step_value=Sim_step.signals.values;
% 
% Sine_bounce     = sine_value(:,1:2);
% Sine_pitch      = sine_value(:,3);
% Step_bounce     = step_value(:,1:2);
% Step_pitch      = step_value(:,3);
% 
% figure
% subplot(2,1,1), plot(Time_sine,rad2deg(Sine_pitch), 'linewidth', 2)
% hold on, xlabel('Time [s]'), ylabel('Displacement [\circ]')
% title('Response of Transfer function at: ' + string(testHZ) + ' Hz')
% subplot(2,1,2), plot(Time_step,rad2deg(Step_pitch), 'linewidth', 2)
% hold on, xlabel('Time [s]'), ylabel('Displacement [\circ]')
% legend('\chi - Pitch')
% 
% sim('Bounce',Time)
% Time_sine=Sim_sine.time; sine_value=Sim_sine.signals.values;
% Time_step=Sim_sine.time; step_value=Sim_step.signals.values;
% 
% Sine_bounce     = sine_value(:,1:2);
% Sine_pitch      = sine_value(:,3);
% Step_bounce     = step_value(:,1:2);
% Step_pitch      = step_value(:,3);
% 
% subplot(2,1,1), plot(Time_sine,rad2deg(Sine_pitch), 'linewidth', 2)
% hold on, xlabel('Time [s]'), ylabel('Displacement [\circ]')
% title('Response of Transfer function at: ' + string(testHZ) + ' Hz')
% subplot(2,1,2), plot(Time_step,rad2deg(Step_pitch), 'linewidth', 2)
% hold on, xlabel('Time [s]'), ylabel('Displacement [\circ]')
% legend('\chi - Passive','\chi - Active')

TF=tf(ss(A,B,C,D));
Tf_bounce=TF(1,:);
Tf_pitch=TF(2,:);
damp(Tf_bounce);
damp(Tf_pitch)

%% Hinf
Time=5
sim('Active_skyhook',Time)
Time_sine=Sim_sine.time; sine_value=Sim_sine.signals.values;
Time_step=Sim_sine.time; step_value=Sim_step.signals.values;

Sine_bounce     = sine_value(:,1:2);
Sine_pitch      = sine_value(:,3);
Step_bounce     = step_value(:,1:2);
Step_pitch      = step_value(:,3);

figure
%
subplot(2,1,1), plot(Time_sine,rad2deg(Sine_pitch), 'linewidth', 2)
hold on, xlabel('Time [s]'), ylabel('Displacement [\circ]')
title('Response of Transfer function at: ' + string(testHZ) + ' Hz')
subplot(2,1,2), plot(Time_step,rad2deg(Step_pitch), 'linewidth', 2)
hold on, xlabel('Time [s]'), ylabel('Displacement [\circ]')
legend('\chi - Pitch')
%

Time_sine=Sim_sine_H.time; sine_value=Sim_sine_H.signals.values;
Time_step=Sim_step_H.time; step_value=Sim_step_H.signals.values;

Sine_bounce     = sine_value(:,1:2);
Sine_pitch      = sine_value(:,3);
Step_bounce     = step_value(:,1:2);
Step_pitch      = step_value(:,3);


subplot(2,1,1), plot(Time_sine,rad2deg(Sine_pitch), 'linewidth', 2)
hold on, xlabel('Time [s]'), ylabel('Displacement [\circ]')
title('Response of Transfer function at: ' + string(testHZ) + ' Hz')
subplot(2,1,2), plot(Time_step,rad2deg(Step_pitch), 'linewidth', 2)
hold on, xlabel('Time [s]'), ylabel('Displacement [\circ]')
legend('\chi - Pitch')
%

sim('Bounce',Time)
Time_sine=Sim_sine.time; sine_value=Sim_sine.signals.values;
Time_step=Sim_sine.time; step_value=Sim_step.signals.values;

Sine_bounce     = sine_value(:,1:2);
Sine_pitch      = sine_value(:,3);
Step_bounce     = step_value(:,1:2);
Step_pitch      = step_value(:,3);

subplot(2,1,1), plot(Time_sine,rad2deg(Sine_pitch), 'linewidth', 2)
hold on, xlabel('Time [s]'), ylabel('Displacement [\circ]')
title('Response of Transfer function at: ' + string(testHZ) + ' Hz')
subplot(2,1,2), plot(Time_step,rad2deg(Step_pitch), 'linewidth', 2)
hold on, xlabel('Time [s]'), ylabel('Displacement [\circ]')
legend('\chi - Passive','\chi - H_{\infty}','\chi - Active')

% %% Bounce
Time = 5
figure
%

sim('Active_skyhook',Time)
Time_sine=Sim_sine.time; sine_value=Sim_sine.signals.values;
Time_step=Sim_sine.time; step_value=Sim_step.signals.values;

Sine_bounce     = sine_value(:,1:2);
Sine_pitch      = sine_value(:,3);
Step_bounce     = step_value(:,1:2);
Step_pitch      = step_value(:,3);

subplot(2,1,1), plot(Time_sine,Sine_bounce, 'linewidth', 2)
hold on, xlabel('Time [s]'), ylabel('Displacement [cm]')
title('Response of Transfer function at: ' + string(testHZ) + ' Hz')
subplot(2,1,2), plot(Time_step,Step_bounce, 'linewidth', 2)
hold on, xlabel('Time [s]'), ylabel('Displacement [cm]')
%

Time_sine=Sim_sine_H.time; sine_value=Sim_sine_H.signals.values;
Time_step=Sim_step_H.time; step_value=Sim_step_H.signals.values;

Sine_bounce     = sine_value(:,2);
Sine_pitch      = sine_value(:,3);
Step_bounce     = step_value(:,2);
Step_pitch      = step_value(:,3);


subplot(2,1,1), plot(Time_sine,Sine_bounce, 'linewidth', 2)
hold on, xlabel('Time [s]'), ylabel('Displacement [cm]')
title('Response of Transfer function at: ' + string(testHZ) + ' Hz')
subplot(2,1,2), plot(Time_step,Step_bounce, 'linewidth', 2)
hold on, xlabel('Time [s]'), ylabel('Displacement [cm]')
%

sim('Bounce',Time)
Time_sine=Sim_sine.time; sine_value=Sim_sine.signals.values;
Time_step=Sim_sine.time; step_value=Sim_step.signals.values;

Sine_bounce     = sine_value(:,2);
Sine_pitch      = sine_value(:,3);
Step_bounce     = step_value(:,2);
Step_pitch      = step_value(:,3);

subplot(2,1,1), plot(Time_sine,Sine_bounce, 'linewidth', 2)
hold on, xlabel('Time [s]'), ylabel('Displacement [cm]')
title('Response of Transfer function at: ' + string(testHZ) + ' Hz')
subplot(2,1,2), plot(Time_step,Step_bounce, 'linewidth', 2)
hold on, xlabel('Time [s]'), ylabel('Displacement [cm]')
legend('impulse','z - Passive','z - H_{\infty}','z - Active')



