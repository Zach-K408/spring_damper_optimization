%% Spring Damper Optimization %% 
function [tt,x,v] = function_eval(m,c,k,f,t,x0,v0,alpha)

% Inputs:
%   m:          Sprung mass [kg] (scalar)
%   c:          Damping coefficient [N*s/m] (scalar)
%   k:          Spring stiffness [N/m] (scalar)
%   f:          Forcing function [N] (anonymous function)
%   t:          Interval of integration [s] (vector)
%   x0:         Initial position [m] (scalar)
%   v0:         Initial velocity [m/s] (scalar)
%   alpha:      Accuracy of ode45() solver
% 
% Outputs:
%   tt:   Output time vector [s] (vector)
%   x:    Position of sprung mass [m] (vector)
%   v:    Velocity of sprung mass [m/s] (vector)

% decomposed 2nd order ODE into a system of 1st order ODEs
z_dot = @(t,z) [z(2);
    (f(t)/m)-(c/m)*z(2)-(k/m)*z(1)];

% solve system of 1st order ODEs
[tt,z] = ode45(z_dot,t,[x0,v0],alpha);

x = z(:,1); % parse the position
v = z(:,2); % parse the velocity

end