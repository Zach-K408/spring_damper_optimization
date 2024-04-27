% pattern search

% define transfer function
transfer_function = @(s,c,k) tf([(c*c2),(k*c2+k2*c),(k*k2)], ...
    [(m1*m2),(m1*c+m1*c2+m2*c),(m1*k+m1*k2+k*m2),(c*c2+c*k2+k*c2),(k*k2)]);

% define bounds of variables
lb = [980, 4300];
ub = [9000, 30000];

% inital point
x0 = [1000,10000];

% Pattern search parameters
stepSize = 10;   % Initial step size
stepFactor = 0.95;   % Step size reduction factor
tolerance = 1e-6;   % Tolerance for convergence
max_iter = 50000;   % Maximum number of iterations

% initialize variables
x = x0;
step = stepSize;
iter = 0;

% start timer
tic;

while step > tolerance && iter < max_iter
    
    % reinitialize arrays
    fs = [];
    candidates = [];
        
    % increment
    iter = iter + 1;
    
    % evaluate objective function at current point
    sys = transfer_function([], x(1), x(2));
    y = lsim(sys,zr,t);
    max_a = max(diff(diff(y)) / (t(2) - t(1))^2);
    f = max_a;
   
    % generate candidate points in the neighborhood
    candidates = [x + step * eye(2); x - step * eye(2)];

    % evaluate feasibility of candidate points
        feasible_candidates = [];
        for i = 1:length(candidates)
            candidate = candidates(i, :);
            
            % check domain constraints
            if all(candidate >= lb) && all(candidate <= ub)
                feasible_candidates = [feasible_candidates; candidate];
            end
        end
    
        % evaluate objective function at candidate points
        for i = 1:length(feasible_candidates)
            sys = transfer_function([], candidates(i,2), candidates(i,1));
            y = lsim(sys, zr, t);
            max_a = max(diff(diff(y)) / (t(2) - t(1))^2);
            fs = [fs,max_a];
        end
     
    % find the index of the best candidate point
    best_index = find(fs == min(fs));
   
    % update the current point
    x = feasible_candidates(best_index,:);
    
    % update the step size
    step = step * stepFactor;
    
end

if iter == max_iter
    disp('Maximum number of iterations reached.');
else

    disp('Converged to the following minimum:');
    disp(x);
    disp('Objective function value at the minimum:');
    disp(f);
end

% Stop time
computation_time3 = toc;

disp(['Computational time: ' num2str(computation_time3) ' seconds']);

