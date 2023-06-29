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

% num2 = [m1*c2,(c1*c2+m1*k2),(k1*c2+c1*k2),k1*k2];
% den2 = [m1*m2,(m1*c1+m1*c2+m2*c1),(m1*k1+m1*k2+k1*m2),(c1*c2+c1*k2+k1*c2),k1*k2];
% 
% sys2 = tf(num2,den2);
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

% f1 = figure('name','Displacement of Sprung Mass');
% lsim(sys1,zr,s/v)
% title('Displacement of Sprung Mass')
% ylabel('Road Displacement (m)')
% legend('Displacement of Sprung Mass')
% ylim([-0.1,0.1])

% f2 = figure('name','Displacement of Unsprung Mass');
% lsim(sys2,zr,s/v)
% ylabel('Road Displacement (m)')
% title('Displacement of Unsprung Mass')

%grid on

%% Generate 3d acceleration plot of sprung mass

% initialize
coordinates = [];
max_responses = [];
max_accelerations = [];

% define transfer function as function handle
transfer_function = @(s,c,k) tf([(c*c2),(k*c2+k2*c),(k*k2)], ...
    [(m1*m2),(m1*c+m1*c2+m2*c),(m1*k+m1*k2+k*m2),(c*c2+c*k2+k*c2),(k*k2)]);

% define domain
c_range = 50:50:5000;
k_range = 1000:500:50000;

% evalute function at each coordinate
for c = c_range
    for k = k_range
        sys = transfer_function([], c, k);
        y = lsim(sys, zr, t);
        max_response = max(abs(y));
        max_acceleration = max(abs(diff(diff(y)) / (t(2) - t(1))^2));
        coordinates = [coordinates; c, k];
        max_responses = [max_responses; max_response];
        max_accelerations = [max_accelerations; max_acceleration];
        
    end
end

% Reshape the coordinates and maximum response values into a grid
[C, K] = meshgrid(c_range, k_range);
Z = reshape(max_accelerations, length(k_range), length(c_range));

% generate 3D plot
f1 = figure('name','Acceleration Method 1');
surf(C,K,Z);
colorbar;
xlabel('Damping Coefficient')
ylabel('Spring Coefficent')
zlabel('Sprung Mass Acceleration')
title('Function to Optimize')

% create a table coordinates
table = array2table([coordinates, max_responses], "VariableNames",{'c','k','z'});

disp(table)

%% Generate 3d acceleration plot of sprung mass (method 2)

% initialize
coordinates = [];
max_responses = [];
max_accelerations = [];

% define transfer function as function handle
transfer_function1 = @(s,c,k) tf([(c*c2),(k*c2+k2*c),(k*k2)], ...
    [(m1*m2),(m1*c+m1*c2+m2*c),(m1*k+m1*k2+k*m2),(c*c2+c*k2+k*c2),(k*k2)]);

transfer_function2 = @(s,c,k) tf([m1*c2,(c*c2+m1*k2),(k*c2+c*k2),k*k2], ...
    [m1*m2,(m1*c+m1*c2+m2*c),(m1*k+m1*k2+k*m2),(c*c2+c*k2+k*c2),k*k2]);

% define domain
c_range = 50:50:5000;
k_range = 1000:500:50000;

% evalute function at each coordinate
for c = c_range
    for k = k_range
        sys1 = transfer_function1([], c, k);
        y1 = lsim(sys1, zr, t);
        sys2 = transfer_function2([], c, k);
        y2 = lsim(sys2, zr, t);
        z = y1(2:end)-y2(2:end);
        z1_dot = diff(y1) / (t(2)-t(1));
        z2_dot = diff(y1) / (t(2)-t(1));

        a = (k*z + c*(z1_dot-z2_dot))/m1;
        
        max_response = max(abs(z));

        max_acceleration = max(abs(a));
        coordinates = [coordinates; c, k];
        max_responses = [max_responses; max_response];
        max_accelerations = [max_accelerations; max_acceleration];
        
    end
end

% Reshape the coordinates and maximum response values into a grid
[C, K] = meshgrid(c_range, k_range);
Z = reshape(max_accelerations, length(k_range), length(c_range));

% generate 3D plot
f2 = figure('name','Acceleration Method 2');
surf(C,K,Z);
colorbar;
xlabel('Damping Coefficient')
ylabel('Spring Coefficent')
zlabel('Sprung Mass Acceleration')
title('Function to Optimize')

% create a table coordinates
table = array2table([coordinates, max_responses], "VariableNames",{'c','k','z'});

disp(table)
