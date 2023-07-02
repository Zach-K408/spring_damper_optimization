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

m1 = 500; % [kg]
m2 = 30; % [kg]
k1 = 20000; % [N/m]
k2 = 200000; % [N/m]
c1 = 1500; % [Ns/m]
c2 = 1; % [Ns/m]

%% Laplacian model
% 
% num = [(c1*c2),(k1*c2+k2*c1),(k1*k2)];
% den = [(m1*m2),(m1*c1+m1*c2+m2*c1),(m1*k1+m1*k2+k1*m2),(c1*c2+c1*k2+k1*c2),(k1*k2)];
% 
% sys1 = tf(num, den);
% display(sys1)
%
% num2 = [m1*c2,(c1*c2+m1*k2),(k1*c2+c1*k2),k1*k2];
% den2 = [m1*m2,(m1*c1+m1*c2+m2*c1),(m1*k1+m1*k2+k1*m2),(c1*c2+c1*k2+k1*c2),k1*k2];
% 
% sys2 = tf(num2,den2);
%display(sys2)

%% Road Profile

N = 1000; % Number of samples

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

%
f0 = figure('name','Road Profile');
plot(t,zr)
xlabel('Time (s)')
ylabel('Road Displacement (m)')
title('Road Profile')

%% Generate 3d acceleration plot of sprung mass

% initialize
coordinates = [];
max_displacements = [];
max_accelerations = [];

% define transfer function as function handle
transfer_function = @(s,c,k) tf([(c*c2),(k*c2+k2*c),(k*k2)], ...
    [(m1*m2),(m1*c+m1*c2+m2*c),(m1*k+m1*k2+k*m2),(c*c2+c*k2+k*c2),(k*k2)]);

% define domain
c_range = 980:25:4300;
k_range = 9000:100:30000;

% evalute function at each coordinate
for c = c_range
    for k = k_range
        sys = transfer_function([], c, k);
        y = lsim(sys, zr, t);
        max_displacement = max(y);
        max_acceleration = max(diff(diff(y)) / (t(2) - t(1))^2);
        coordinates = [coordinates; c, k];
        max_displacements = [max_displacements; max_displacement];
        max_accelerations = [max_accelerations; max_acceleration];
        
    end
end

% Reshape the coordinates and maximum response values into a grid
[C, K] = meshgrid(c_range, k_range);
Z = reshape(max_accelerations, length(k_range), length(c_range));

% generate 3D plot
f1 = figure('name','Objective Function');
surf(C,K,Z);
colorbar;
xlabel('Damping Coefficient (Ns/m)')
ylabel('Spring Coefficent (N/m)')
zlabel('Sprung Mass Acceleration (m/s^2)')
title('Surface Plot of Objective Function')

%%
% generate 2D contour plot
f2 = figure('name','Objective Function');
contourf(C,K,Z,100);
colormap();
colorbar();
xlabel('Damping Coefficient (Ns/m)')
ylabel('Spring Coefficent (N/m)')
zlabel('Sprung Mass Acceleration (m/s^2)')
title('Contour Plot of Objective Function')

%%

% create a table coordinates
table = array2table([coordinates, max_accelerations], "VariableNames",{'c','k','accel'});

writetable(table,'output.csv');
%% Generate 3d acceleration plot of sprung mass (method 2)
% 
% % initialize
% coordinates = [];
% max_relative_displacements = [];
% max_accelerations = [];
% 
% % define transfer function as function handle
% transfer_function1 = @(s,c,k) tf([(c*c2),(k*c2+k2*c),(k*k2)], ...
%     [(m1*m2),(m1*c+m1*c2+m2*c),(m1*k+m1*k2+k*m2),(c*c2+c*k2+k*c2),(k*k2)]);
% 
% transfer_function2 = @(s,c,k) tf([m1*c2,(c*c2+m1*k2),(k*c2+c*k2),k*k2], ...
%     [m1*m2,(m1*c+m1*c2+m2*c),(m1*k+m1*k2+k*m2),(c*c2+c*k2+k*c2),k*k2]);
% 
% % define domain
% c_range = 980:25:4300;
% k_range = 9000:100:30000;
% 
% % evalute function at each coordinate
% for c = c_range
%     for k = k_range
%         sys1 = transfer_function1([], c, k);
%         y1 = lsim(sys1, zr, t);
%         sys2 = transfer_function2([], c, k);
%         y2 = lsim(sys2, zr, t);       
% 
%         z = y1-y2;
% 
%         max_relative_displacement = max(z);
%         relative_velocity = diff(z)/(t(2)-t(1));
%        
%         max_acceleration = max((k*z(2:end)+c*relative_velocity)/m1);
% 
%         coordinates = [coordinates; c, k];
%         max_relative_displacements = [max_relative_displacements; max_relative_displacement];
%         max_accelerations = [max_accelerations; max_acceleration];
% 
% 
% 
%     end
% end
% 
% 
% % Reshape the coordinates and maximum response values into a grid
% [C, K] = meshgrid(c_range, k_range);
% Z = reshape(max_accelerations, length(k_range), length(c_range));
% 
% % generate 3D plot
% f2 = figure('name','Acceleration Method 2');
% surf(C,K,Z);
% colorbar;
% xlabel('Damping Coefficient')
% ylabel('Spring Coefficent')
% zlabel('Sprung Mass Acceleration')
% title('Function to Optimize')
% 
% % create a table coordinates
% table = array2table([coordinates, max_relative_displacements], "VariableNames",{'c','k','z**'});
% 
% disp(table)
