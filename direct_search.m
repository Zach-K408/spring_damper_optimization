% Brute Force


% add a feature to record the computational time

data = readtable('output.csv');

arr = table2array(data);

X = arr(:, 1);
Y = arr(:, 2);
Z = arr(:, 3);

minZ = Z(1);  
minIndex = 1;  

for i = 2:numel(Z)
    if Z(i) < minZ
        minZ = Z(i);
        minIndex = i;
    end
end

minX = X(minIndex);
minY = Y(minIndex);

disp('The minimum is found at the coordinates:');
disp([minX, minY, minZ]);
