% Nelder-Mead Method

% objective function Z(X, Y)
transfer_function = @(s,c,k) tf([(c*c2),(k*c2+k2*c),(k*k2)], ...
    [(m1*m2),(m1*c+m1*c2+m2*c),(m1*k+m1*k2+k*m2),(c*c2+c*k2+k*c2),(k*k2)]);

N = 2; % number of variables
alpha = 1; % reflection coefficient
beta = 0.5; % contraction coefficient
gamma = 2; % expansion coefficient
tol = 1e-6; % termination tolerance
max_iter = 100; % maximum number of iterations

% test domain for X and Y
c_min = 980;
c_max = 4300;
k_min = 9000;
k_max = 30000;

% initialize simplex
x0 = [1000, 10000]; % initial point (can be any valid point within the test domain)
simplex = zeros(N+1, N);
simplex(1,:) = x0;

for i = 1:N
    point = x0;
    point(i) = point(i) + 1;
    simplex(i+1,:) = point;
end

% start timer
tic;

% Nelder-Mead iterations
iter = 0;
while iter < max_iter
    % evaluate objective function at simplex points
    values = zeros(N+1, 1);
    for i = 1:N+1
        sys = transfer_function([], simplex(i,1), simplex(i,2));
        y = lsim(sys,zr,t);
        max_a = max(diff(diff(y)) / (t(2) - t(1))^2);
        values(i) = max_a;
    end
    
    % sort simplex based on objective function values
    [values, order] = sort(values);
    simplex = simplex(order,:);
    
    % check termination condition
    if max(abs(values(2:end) - values(1))) < tol
        break;
    end
    
    % calculate centroid (excluding worst point)
    centroid = mean(simplex(1:N,:));
    
    % reflection
    xr = centroid + alpha * (centroid - simplex(end,:)); % alpha reflection coordinates
    xr = max(min(xr, [c_max, k_max]), [c_min, k_min]); % apply domain constraints
    % function evaluation at alpha reflection
    sys = transfer_function([], xr(1), xr(2));
    y = lsim(sys,zr,t);
    max_a = max(diff(diff(y)) / (t(2) - t(1))^2);
    fr = max_a;
    
    if fr < values(1) % if alpha relection value < lowest point of previous simplex
        % Expansion
        xe = centroid + gamma * (xr - centroid); % gamma reflection coordinates
        xe = max(min(xe, [c_max, k_max]), [c_min, k_min]); % apply domain constraints
        % function evaluation at gamma reflection
        sys = transfer_function([], xe(1), xe(2));
        y = lsim(sys,zr,t);
        max_a = max(diff(diff(y)) / (t(2) - t(1))^2);
        fe = max_a; 
        
        if fe < fr % if gamma reflection less than alpha reflection
            simplex(end,:) = xe; % replace with gamma reflection
        else
            simplex(end,:) = xr; % otherwise replace with alpha reflection
        end
    else % if alpha reflection > lowest point of previous simplex 
        if fr < values(end-1) % if alpha reflection is < second lowest point
            simplex(end,:) = xr; % replace with alpha reflection
        else
            % Contraction
            if fr < values(end) % if alpha reflection < its previous point
                xc = centroid + beta * (xr - centroid); % beta reflection coordinate
                xc = max(min(xc, [c_max, k_max]), [c_min, k_min]); % apply domain constraints
                % function evaluation at beta reflection
                sys = transfer_function([], xc(1), xc(2));
                y = lsim(sys,zr,t);
                max_a = max(diff(diff(y)) / (t(2) - t(1))^2);
                fc = max_a;
                
                if fc < fr % if beta reflection < alhpa reflection
                    simplex(end,:) = xc; % replace with beta reflection
                else % if beta reflection is worse than alpha reflection
                    % shrink
                    for i = 2:N+1
                        simplex(i,:) = simplex(1,:) + 0.5 * (simplex(i,:) - simplex(1,:)); % bring the other two points 0.5 times closer to the best point
                        simplex(i,:) = max(min(simplex(i,:), [c_max, k_max]), [c_min, k_min]); % apply domain constraints
                    end
                end
            else
                % shrink
                for i = 2:N+1
                    simplex(i,:) = simplex(1,:) + 0.5 * (simplex(i,:) - simplex(1,:)); % bring the other two points 0.5 times closer to the best point
                    simplex(i,:) = max(min(simplex(i,:), [c_max, k_max]), [c_min, k_min]); % apply domain constraints
                end
            end
        end
    end
    
    iter = iter + 1;
end

% stop time
computation_time2 = toc;

% extract minimum point
X = simplex(1,1);
Y = simplex(1,2);
sys = transfer_function([], simplex(1,1), simplex(1,2));
y = lsim(sys,zr,t);
max_a = max(diff(diff(y)) / (t(2) - t(1))^2);
min_Z = max_a;

% display results
disp(['Minimum found at (X, Y) = (' num2str(X) ', ' num2str(Y) ')']);
disp(['Minimum value of Z: ' num2str(min_Z)]);
disp(['Computational time: ' num2str(computation_time2) ' seconds']);
