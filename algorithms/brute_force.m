% Direct Search Method

% initialize
coordinates = [];
max_accelerations = [];

% define transfer function as function handle
transfer_function = @(s,c,k) tf([(c*c2),(k*c2+k2*c),(k*k2)], ...
    [(m1*m2),(m1*c+m1*c2+m2*c),(m1*k+m1*k2+k*m2),(c*c2+c*k2+k*c2),(k*k2)]);

% define domain
c_range = 980:25:4300;
k_range = 9000:100:30000;

% start timer
tic;

% evalute function at each coordinate
for c = c_range
    for k = k_range
        sys = transfer_function([], c, k);
        y = lsim(sys, zr, t);
        max_acceleration = max(diff(diff(y)) / (t(2) - t(1))^2);
        coordinates = [coordinates; c, k];
        max_accelerations = [max_accelerations; max_acceleration];
        
    end
end

% create a table coordinates
table = array2table([coordinates, max_accelerations], "VariableNames",{'X','Y','Z'});

% convert table into array for indexing
data = table2array(table);

X = data(:,1);
Y = data(:,2);
Z = data(:,3);

minZ = Z(1);  
minIndex = 1;  

% search for minimum Z value in array and save the index
for i = 2:numel(Z)
    if Z(i) < minZ
        minZ = Z(i);
        minIndex = i;
    end
end

% stop timer
computational_time1 = toc;

minX = X(minIndex);
minY = Y(minIndex);

disp('The minimum is found at the coordinates:');
format long
disp([minX, minY, minZ]);
disp(['Computational time: ' num2str(computational_time1) ' seconds']);