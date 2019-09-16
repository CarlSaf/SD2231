%% This is a Matlab file for designing H_infinity controller (assignment 3
%% of SD2231)
% clear all
s=tf('s');

%Please run bounce_controll first!!

% systme parameters
m=22000;   %kg
j=700e3;   %kgm^2
c=40e3;    %Ns/m
k=2*300e3*(1-0.15); %N/m
L=6;       %m
% testHZ=1
% sinfreq=testHZ*2*pi;

%% State space model for skyhook contorl
Ask=[0 1 0 0
    -2*k/m 0 0 0
    0 0 0 1
    0 0 -2*k*L^2/j 0];
Bsk=[0 0 0 0
    k/m k/m -1/m -1/m
    0 0 0 0
    -L*k/j L*k/j L/j -L/j];
Csk=[0 1 0 0
    0 0 0 1];
Dsk=zeros(2,4);

%% H_inf using linmod syntax

%state space: The same as skyhook

      
%Weighting functions

%For penalizing actuator force
Wa1=(0.00175*s+1)/(0.00025*s+1);
Wa2=Wa1;

%For penalizing bounce and pitch motions
eps=1;
wnb=7.86;            %Find the right equation or value for wnb
wnchi=7.86;          %Find the right equation or value for wnchi
s1b=-eps+1i*sqrt(wnb^2-eps^2);
s2b=-eps-1i*sqrt(wnb^2-eps^2);
s1chi=-eps+1i*sqrt(wnchi^2-eps^2);
s2chi=-eps-1i*sqrt(wnchi^2-eps^2);
kb=5.5e3; 
kchi=2e4;
Wb=(kb*s1b*s2b)/((s-s1b)*(s-s2b));
Wchi=(kchi*s1chi*s2chi)/((s-s1chi)*(s-s2chi));

%Extracting the extended model
[A_Pe,B_Pe,C_Pe,D_Pe] = linmod('Extended_model');% state space parameters of the extended system: Pe
Pe=ss(A_Pe,B_Pe,C_Pe,D_Pe);

%Calculating the controller
ncont = 2;%Number of control inputs
nmeas = 2;%Number of measured outputs provided to the controller
Pe=minreal(Pe);%This syntax cancels pole-zero pairs in transfer
%functions. The output system has minimal order and the same response
%characteristics as the original model.
[K,Pec,gamma,info]=hinfsyn(Pe,nmeas,ncont,'method','lmi'); % for working with the error
[Ainf, Binf, Cinf, Dinf]=ssdata(K);


sim('H_inf')
%Now use the controller K in your simulation

% 
% figure, bodemag(Wa1, Wa2, Wb, Wchi)
Time_sine=Sim_sine_H.time; sine_value=Sim_sine_H.signals.values;
Time_step=Sim_step_H.time; step_value=Sim_step_H.signals.values;
Sine_bounce     = sine_value(:,2);
Sine_pitch      = sine_value(:,3);
Step_bounce     = step_value(:,2);
Step_pitch      = step_value(:,3);


% subplot(2,1,1), plot(Time_sine,rad2deg(Sine_pitch), 'linewidth', 2)
% hold on, xlabel('Time [s]'), ylabel('Displacement [\circ]')
% title('Response of Transfer function at: ' + string(testHZ) + ' Hz')
% subplot(2,1,2), plot(Time_step,rad2deg(Step_pitch), 'linewidth', 2)
% hold on, xlabel('Time [s]'), ylabel('Displacement [\circ]')
% legend('\chi - Pitch')

subplot(2,1,1), plot(Force.time, Force.signals.values(:,1), 'linewidth', 2)
hold on, xlabel('Time [s]'), ylabel('Actuator Force [N]')
title('Response of Transfer function at: ' + string(testHZ) + ' Hz')
subplot(2,1,2), plot(Time_step,Step_bounce, 'linewidth', 2)
hold on, xlabel('Time [s]'), ylabel('Displacement [cm]')
% legend('+m','-m','+J','-J','+c','-c','+k','-k')

