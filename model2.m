%% Quarter Car Model %% 
clc
clear
close all;

% Inputs:
%   m1:          Sprung mass [kg] (scalar)
%   m2:          Unsprung mass [kg] (scalar)
%   c1:          Damping coefficient of sprung mass [N*s/m] (scalar)
%   c2:          Damping coefficient of unsprung mass [N*s/m] (scalar)
%   k1:          Spring stiffness of sprung mass [N/m] (scalar)
%   k2:          Spring stiffness of unsprung mass [N/m] (scalar)
%   u:          Road profile function [m] (function)
%   t:          Interval of integration [s] (vector)
%   z1:         Displacement of sprung mass [m] (scalar)
%   z1:         Displacement of unsprung mass [m] (scalar)
%   v:         Velocity of vehicle [m/s] (scalar)
%   alpha:      Accuracy of ode45() solver
% 
% Outputs:
%   z1:    Displacement of sprung mass [m] (vector)
%   z2:    Displacement of unsprung mass [m] (vector)

% model params
m1 = 290; % [kg]
m2 = 15; % [kg]
k1 = 16200; % [N/m]
k2 = 191000; % [N/m]
c1 = 1000; % [Ns/m]
c2 = 2500; % [Ns/m]
%t=0:0.001:7; % [s]
%u=0.1; % [m]

num = [(c1*c2),(k1*c2+k2*c1),(k1*k2),0,0];
den = [(m1*m2),(m1*c1+m1*c2+m2*c1),(m1*k1+m1*k2+k1*m2),(c1*c2+c1*k2+k1*c2),(k1*k2)];

sys1 = tf(num, den);
display(sys1)
num2=[m1*c2,(c1*c2+m1*k2),(k1*c2+c1*k2),k1*k2,0,0];
den2=[m1*m2,(m1*c1+m1*c2+m2*c1),(m1*k1+m1*k2+k1*m2),(c1*c2+c1*k2+k1*c2),k1*k2];

sys2 = tf(num2,den2);
display(sys2)

[u,t] = gensig("square",10,20);

f1 = figure('name','Displacement of Sprung Mass');
lsim(sys1,u,t)
title('Displacement of Sprung Mass')
ylabel('Road Displacement (m)')
legend('Displacement of Sprung Mass')
f2 = figure('name','Displacement of Unsprung Mass');
lsim(sys2,u,t)
ylabel('Road Displacement (m)')
title('Displacement of Unsprung Mass')


