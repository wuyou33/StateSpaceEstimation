clc; close all; clearvars; clear memoize;
addpath(genpath('./'));

fprintf('\nThe demonstration of state space esitmation problem of bearing and frequency tracking problem.\n');

% {'srukf', 'pf', 'gspf', 'sppf', 'gmsppf'}
filter_type = 'gmsppf';

model = gssm_bearing_and_frequency_tracking('init');

arg.model = model;
arg.type = 'state';
arg.tag = 'State estimation for bearings and frequency tracking problem';
inference_model = inference_model_generator(arg);

% max. time k=1..T
T = 100;
state_noise = model.processNoise.sample(model.processNoise, T);
observ_noise = model.observationNoise.sample(model.observationNoise, T);

%Generate the observer (submarine) trajectory

%Initial state of submarine
sub = zeros(inference_model.stateDimension, T);

sub(1, 1) = 0;
sub(2, 1) = 10;
sub(3, 1) = 0;
sub(4, 1) = 0;
sub(5, 1) = 350;

for k = 2:(T/2-1);
    sub(:,k) = model.transition_fun(model, sub(:, k-1), state_noise(:, k-1), []);
end

% At mid-course, the sub changes its course (90 grad turn)
sub_speedX = sub(2, T/2-1);
sub_speedY = sub(4, T/2-1);
sub(2, T/2-1) = sub_speedY;
sub(4, T/2-1) = sub_speedX;

for k = T/2:T;
    sub(:,k) = model.transition_fun(model, sub(:, k-1), state_noise(:, k-1), []);
end

% Generate the target trajectory
x_true = zeros(inference_model.stateDimension, T);
z = zeros(inference_model.observationDimension, T);


range_0        = 2000 + 100*randn(1);
bearing_0      = -pi + 2*pi*rand(1);
frequency_0    = 300;
course_0       = -pi + 2*pi*rand(1);
speed_0        = 12+1*randn(1);

x_true(:, 1)=[range_0.*cos(bearing_0) + sub(1,1);
    speed_0.*cos(course_0);
    range_0.*sin(bearing_0) + sub(1,1);
    speed_0.*sin(course_0);
    frequency_0];

z(:, 1) = model.observation_fun( model, x_true(:,1), observ_noise(:,1), sub(:,1));
for k = 2:T
    x_true(:, k) = model.transition_fun(model, x_true(:, k-1), state_noise(:,k-1), []);
    z(:, k) = model.observation_fun(model, x_true(:, k), observ_noise(:,k), sub(:,k));
end

true_range   = sqrt((x_true(1, :) - sub(1, :)).^2 + (x_true(3, :) - sub(3, :)).^2);
true_bearing = atan2(x_true(3, :) - sub(3, :), x_true(1, :) - sub(1, :));
true_frequency = zeros(1, T);
for m = 1:T
    true_frequency(1,m) = x_true(5, m) * (1+1/1500*((sub(2, m) - x_true(2, m)) * cos(true_bearing(1, m))+(sub(4, m)-x_true(4, m)) * sin(true_bearing(1, m))));
end

%--- Display target trajectory detail

figure(1);
clf;
p1 = plot(x_true(1, :), x_true(3, :),'r-');
hold on;
p2 = plot(x_true(1, 1), x_true(3, 1),'g*');
p3 = plot(x_true(1, end), x_true(3, end),'k*');

p4 = plot(sub(1,:),sub(3,:),'b-');
hold on;
p5 = plot(sub(1,1), sub(3,1), 'g*');
p6 = plot(sub(1, end), sub(3, end),'k*');
hold off;
legend([p1 p2 p3 p4],'target trajectory','position : k=0', ['position : k=' num2str(T)], 'trajectory');
xlabel('x');
ylabel('y');
title('Target and submarine trajectory');
grid on;

figure(2);
subplot(211);
p11 = plot(1 : T, true_range, 'b-o');
xlabel('k');
ylabel('Distance');
title('Distance');
legend(p11,'true');
subplot(212);
p13 = plot(1:T,true_bearing,'b-o');
hold on;
p14 = plot(1:T,z(1,:),'k+');
hold off;
legend([p13 p14], 'true bearing','measured bearing');
xlabel('time : k');
ylabel('bearing : radians');
title('Bearing');

figure(3);
p17 = plot(true_frequency(1, :), 'b-o');
hold on;
p18 = plot(z(2, :), 'k+');
hold off;
legend([p17 p18], 'true frequency', 'measured frequency');
xlabel('time : k');
ylabel('frequency : Hertz');
title('Doppler shifted frequency');
axis([0 T 295 305]);
drawnow;

% Setup estimation buffers
x_est = zeros(inference_model.stateDimension, T);
x_cov_est = eye(inference_model.stateDimension);

% Determine initial uncertainty in vehicle position
Nstat = 10000;

bearing_stat      = true_bearing(1) + sqrt(model.observationNoise.covariance(1, 1)) * randn(1, Nstat);
range_stat        = 2000+500*randn(1, Nstat);
course_stat       = course_0 + 2*rand(1, Nstat);
speed_stat        = 12 + randn(1, Nstat);
frequency_stat    = 300 + sqrt(model.observationNoise.covariance(2, 2)) * randn(1, Nstat);

Xstat = [range_stat.*cos(bearing_stat) + sub(1, 1);
    speed_stat.*cos(course_stat);
    range_stat.*sin(bearing_stat) + sub(1, 1);
    speed_stat.*sin(course_stat);
    frequency_stat];

Mu0 = mean(Xstat, 2);
P0  = cov(Xstat');

switch filter_type
    case 'pf'
        numParticles = 5e5;
    case 'gspf'
        numParticles = 5e4;
    case 'gmsppf'
        numParticles = 5e3;
    case 'sppf'
        numParticles = 7e2;
    otherwise
        numParticles = 1e3;
end

bearing_stat      = true_bearing(1) + sqrt(model.observationNoise.covariance(1, 1)) * randn(1, numParticles);
range_stat        = 2000+500*randn(1, numParticles);
course_stat       = course_0 + 2 * rand(1, numParticles);
speed_stat        = 12 + randn(1, numParticles);
frequency_stat    = 300 + sqrt(model.observationNoise.covariance(2, 2)) * randn(1, numParticles);

initialParticles = [range_stat.*cos(bearing_stat) + sub(1,1);
    speed_stat.*cos(course_stat);
    range_stat.*sin(bearing_stat) + sub(1,1);
    speed_stat.*sin(course_stat);
    frequency_stat];

[proc_noise_infer, obs_noise_infer, inference_model] = inference_noise_generator(inference_model, filter_type);

tic;
switch filter_type
    case 'pf'
        particle_set.particlesNum = numParticles;
        particle_set.particles = initialParticles;
        particle_set.weights = (1/numParticles)*ones(1, numParticles);
        
        inference_model.resampleThreshold = 1;
        inference_model.estimateType = 'mean';
        inference_model.resampleMethod = 'systematic';
        for k = 2 : T
            [x_est(:, k), particle_set] = pf(particle_set, proc_noise_infer, obs_noise_infer, z(:, k), [], sub(:, k), inference_model);
        end
    case 'gspf'
        particle_set.particlesNum = numParticles;
        particle_set.stateGMM = gmm_fit(initialParticles, 7, [0.001 10], 'sqrt', 1);
        
        inference_model.estimateType = 'mean';
        inference_model.resampleMethod = 'systematic';
        for k = 2 : T
            [x_est(:, k), particle_set] = gspf(particle_set, proc_noise_infer, obs_noise_infer, z(:, k), [], sub(:, k), inference_model);
        end
    case 'gmsppf'
        particle_set.particlesNum = numParticles;
        particle_set.stateGMM = gmm_fit(initialParticles, 7, [0.001 10], 'sqrt', 1);
        inference_model.estimateType = 'mean';
        inference_model.resampleMethod = 'systematic';
        % inference_model.spkfType = 'srukf';
        % inference_model.spkfParams  = [1 2 0];
        inference_model.spkfType = 'srcdkf';
        inference_model.spkfParams = sqrt(3);
        
        for k = 2 : T
            [x_est(:, k), particle_set] = gmsppf(particle_set, proc_noise_infer, obs_noise_infer, z(:, k), [], sub(:, k), inference_model);
        end
    case 'sppf'
        inference_model.estimateType = 'mean';    % estimate type for Xh
        inference_model.resampleMethod = 'residual';
        % inference_model.spkfType = 'srukf';
        % inference_model.spkfParams  = [1 2 0];
        inference_model.spkfType = 'srcdkf';
        inference_model.spkfParams = sqrt(3);
        inference_model.resampleThreshold = 1;
        
        % generate Gaussian system noise sources for internal SPKFs
        [proc_noise_gauss, observ_noise_gauss, ~] = inference_noise_generator(inference_model, inference_model.spkfType);
        
        % build ParticleFilter data structure
        particle_set.particlesNum = numParticles;
        particle_set.particles = initialParticles;
        particle_set.particlesCov = repmat(chol(P0)', [1 1 numParticles]);
        particle_set.processNoise = proc_noise_gauss;
        particle_set.observationNoise = observ_noise_gauss;
        particle_set.weights = column_vector_replicate(1 / numParticles,numParticles);
        
        for k = 2 : T
            [x_est(:, k), particle_set] = sppf(particle_set, proc_noise_infer, obs_noise_infer, z(:, k), [], sub(:, k), inference_model);
        end
    case 'srukf'
        inference_model.spkfParams  = [1 2 0];
        x_cov_est = chol(P0)';
        x_est(:, 1) = Mu0;
        
        for k = 2 : T
            [x_est(:, k), x_cov_est] = srukf(x_est(:, k-1), x_cov_est, proc_noise_infer, obs_noise_infer, z(:, k), inference_model, [], sub(:, k));
        end
    otherwise
        error([' Algorithm type is not supported''' filter_type '''']);
end
fprintf('%s done', filter_type);
toc;
fprintf('\n\n');


range_estimate = sqrt((x_est(1, :) - sub(1, :)).^2 + (x_est(3, :)-sub(3, :)).^2);
bearing_estimate = atan2(x_est(3, :) - sub(3, :), x_est(1, :) - sub(1, :));
frequency_estimate = zeros(1, T);
for m = 1:T
    frequency_estimate(1, m) = x_est(5, m)*(1+1/1500*((sub(2, m)-x_est(2, m))*cos(bearing_estimate(1, m))+(sub(4, m)-x_est(4, m))*sin(bearing_estimate(1, m))));
end

range_error     =  range_estimate - true_range;
bearing_error   =  bearing_estimate - true_bearing;
pos_error       =  sqrt((x_est([1; 3], :) - x_true([1; 3], :)).^2);


figure(1);
hold on;
p7=plot(x_est(1, :), x_est(3, :), 'r.', 'LineWidth', 1);
legend([p1 p2 p3 p4 p7], 'true target trajectory','position : k=0', ['position : k=' num2str(T)], 'submarine trajectory','estimated target''s trajectory');
xlabel('x');
ylabel('y');
title('Target and submarine trajectory');
hold off;

figure(2);
subplot(211);
hold on;
p12 = plot(1:T, range_estimate, 'r-', 'LineWidth', 1);
xlabel('k');
ylabel('range');
title('Range Profile');
grid on;
legend([p11 p12], 'true', 'inferred');
hold off;
subplot(212);
hold on;
p15 = plot(1:T, bearing_estimate, 'r-', 'LineWidth', 1);
xlabel('t');
ylabel('bearing');
title('Bearing Profile');
legend([p13 p14 p15], 'true', 'measured', 'inferred');
hold off;
grid on;

figure(3);
p16 = plot(x_est(5, :), 'r-', 'LineWidth', 1);
hold on;
p17 = plot(x_true(5, :), 'b-', 'LineWidth', 1);
xlabel('t');
ylabel('Source frequency : Hertz');
title('Source frequency profile');
legend([p16 p17], 'estimated frequency', 'true source frequency');
axis([0 T 280 320]);
hold off;
grid on;

figure(4);
p18 = plot(true_frequency(1, :), 'b-o', 'LineWidth', 1);
hold on;
p19 = plot(frequency_estimate(1,:), 'r-', 'LineWidth', 1);
p20 = plot(z(2,:), 'k+', 'LineWidth', 1);
legend([p18 p19 p20], 'true doppler shifted frequency', 'estimated doppler shifted frequency', 'measured doppler shifted frequency');
xlabel('t');
ylabel('frequency : Hertz');
title('Doppler shifted frequency profile');
axis([0 T 295 305]);
hold off;
grid on;

figure(5);
subplot(211);
p21 = plot(atan2(x_true(4, :), x_true(2, :)), 'b-', 'LineWidth', 1);
hold on;
grid on;
p22 = plot(atan2(x_est(4, :), x_est(2, :)), 'r-', 'LineWidth', 1);
hold off;
grid on;
xlabel('t');
ylabel('Course : rad');
title('Course profile (true and estimated)');
legend([p21 p22],'True course','Estimated course');

subplot(212);
p23 = plot(sqrt(x_true(4, :).^2 + x_true(2, :).^2), 'b-', 'LineWidth', 1);
hold on;
grid on;
p24 = plot(sqrt(x_est(4, :).^2 + x_est(2, :).^2), 'r-', 'LineWidth', 1);
hold off;
grid on;
xlabel('t');
ylabel('Speed : m/s');
title('Speed profile (true and estimated)');
legend([p23 p24],'True speed','Estimated speed');


figure(6);
subplot(211);
plot(abs( atan2(x_true(4, :), x_true(2, :)) - atan2(x_est(4, :), x_est(2, :)) ), 'LineWidth', 1);
grid on;
xlabel('t');
ylabel('Course estimation error: rad');
title('Course profile error');

subplot(212);
plot(abs( (sqrt(x_true(4, :).^2 + x_true(2, :).^2)) - (sqrt(x_est(4, :).^2 + x_est(2, :).^2) ) ), 'LineWidth', 1);
grid on;
xlabel('t');
ylabel('Speed estimation error: m/s');
title('Speed profile error');
