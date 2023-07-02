% Genetic Method

% Objective function Z(X, Y)
transfer_function = @(s,c,k) tf([(c*c2),(k*c2+k2*c),(k*k2)], ...
    [(m1*m2),(m1*c+m1*c2+m2*c),(m1*k+m1*k2+k*m2),(c*c2+c*k2+k*c2),(k*k2)]);

% test domain for c and k
lb = [980,9000];
ub = [4300,30000];

% Run the Genetic Algorithm
[x_opt, fval] = ga(@objective, 2, [], [], [], [], lb, ub);

% Display Results
disp('Optimization Results:');
disp(['x_opt = ', num2str(x_opt)]);
disp(['fval = ', num2str(fval)]);

function f = objective(c,k)

    sys = transfer_function([],c,k);
    y = lsim(sys,zr,t);
    f = max(diff(diff(y)) / (t(2) - t(1))^2);
    
end
