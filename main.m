%% Spring Damper Optimization %%
clc
clear
close all;

% parameters
m = 100; % kg
c = 2; % kg/s
k = 15; % kg/s^2

%%write a function that computes the acceleration of the sprung mass when varying the decision
%%variables, namely the stiffness of the spring and damper.

function a = analysis(m,c,k,f)

a = 1/m * (f - c*v - k*x)

end