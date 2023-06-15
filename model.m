%% Quarter Car Model %% 
clc
clear
close all;
%% Model Parameters
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
%   v:          Velocity of vehicle [m/s] (scalar)
%   alpha:      Accuracy of ode45() solver
% 
% Outputs:
%   z1:    Displacement of sprung mass [m] (vector)
%   z2:    Displacement of unsprung mass [m] (vector)

m1 = 290; % [kg]
m2 = 25; % [kg]
k1 = 16200; % [N/m]
k2 = 191000; % [N/m]
c1 = 1000; % [Ns/m]
c2 = 2500; % [Ns/m]

%% Laplacian model

num = [(c1*c2),(k1*c2+k2*c1),(k1*k2)];
den = [(m1*m2),(m1*c1+m1*c2+m2*c1),(m1*k1+m1*k2+k1*m2),(c1*c2+c1*k2+k1*c2),(k1*k2)];

sys1 = tf(num, den);
display(sys1)

num2 = [m1*c2,(c1*c2+m1*k2),(k1*c2+c1*k2),k1*k2];
den2 = [m1*m2,(m1*c1+m1*c2+m2*c1),(m1*k1+m1*k2+k1*m2),(c1*c2+c1*k2+k1*c2),k1*k2];

sys2 = tf(num2,den2);
%display(sys2)

%% Road Profile

N = 500; % Number of samples

v = 80/3.6; % Speed
t = linspace(0, 250/(v), N); % Time vector

s = v*t; % Space coordinate (m)

Om_min = 2*pi/100; % Min frequency
Om_max =2*pi*10; % Max frequency

dOm = (Om_max-Om_min)/(N-1); % Space between frequency samples
Om = Om_min:dOm:Om_max; % Frequency vector
Om_0 = 1; % Ref wavenumber

w = 2; % Waviness
Phi_0 = 2*10e-6; % Depends on the class of the road

Phi = Phi_0.*(Om./Om_0).^(-w); 

rng("default");
Psi = 2*pi*rand(size(Om)); % Phase angles

Amps = sqrt(2*Phi*dOm); % Amplitudes

% Giving the parameters into a structure
p.Amp = Amps;
p.v = v;
p.Om = Om;
p.Psi = Psi;

zr = zeros(size(t)); % Road elevation vector

% Compute the road elevation for each time step

for i=1:length(t)
    zr(i) = road_profile(t(i), p);
end

%% Visualization

f1 = figure('name','Displacement of Sprung Mass');
lsim(sys1,zr,s/v)
title('Displacement of Sprung Mass')
ylabel('Road Displacement (m)')
legend('Displacement of Sprung Mass')
ylim([-0.1,0.1])

% f2 = figure('name','Displacement of Unsprung Mass');
% lsim(sys2,zr,s/v)
% ylabel('Road Displacement (m)')
% title('Displacement of Unsprung Mass')

grid on

