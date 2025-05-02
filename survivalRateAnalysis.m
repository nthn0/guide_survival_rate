%survivalRateAnalysis.m
tic;

%setup
p = gcp('nocreate');
if isempty(p)
    p = parpool('local', min(5, feature('numcores')));
end

[r_dif_1, v_r_1, y_1, v_y_1] = get_input_coords();
n_total_particles = numel(r_dif_1);

r_dif = gpuArray(r_dif_1);
v_r = gpuArray(v_r_1);
y = gpuArray(y_1);
v_y = gpuArray(v_y_1);


%calculations
r_tube = 0.00381;    %m, 0.15 in, radius of inside of tube
m = 1.178983541e-25; %kg, 71 amu
dt = 1e-6;           %s, time step

v_longs = 90:1.875:120;      % 16 values
R_tubes = 0.48:0.005:0.54;   % 12 values
[v_grid, R_grid] = meshgrid(v_longs, R_tubes);
num_points = numel(v_grid);
survival_rate_1d = zeros(numel(v_grid), 1);

parfor i = 1:num_points
    survival_rate_1d(i) = get_survival_rate(r_dif, v_r, y, v_y, v_grid(i), R_grid(i), r_tube, m, dt);
end
survival_rate_grid = reshape(survival_rate_1d, size(v_grid));


%plot
figure(1);
clf;

imagesc(v_longs, R_tubes, survival_rate_grid);
set(gca, 'YDir', 'normal');
colormap('hot');
c = colorbar;
c.Label.String = 'Survival Rate (%)';
title('Survival Rate for Different Radii and Initial Velocities, 71 amu');
subtitle('Computed with dt = ' + string(dt) + ', ' + string(n_total_particles) + ' particles' );
xlabel('Initial Longitudinal Velocity (m/s)'); 
ylabel('Guide Radius (m)');


%printout - can copy and paste into Python
fprintf('v_longs = [%.6f', v_grid(1));
for i = 2:num_points
    fprintf(', %.6f', v_grid(i));
end
fprintf(']\nR_tubes = [%.6f', R_grid(1));
for i = 2:num_points
    fprintf(', %.6f', R_grid(i));
end
fprintf(']\nsurvival_rates = [%.2f', survival_rate_1d(1));
for i = 2:num_points
    fprintf(', %.2f', survival_rate_1d(i));
end
fprintf(']\n');

t = toc;
fprintf('Elapsed time %s\n', string(duration(0, 0, t, 'Format', 'hh:mm:ss')));
