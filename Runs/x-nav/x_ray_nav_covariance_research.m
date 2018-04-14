close all; clc; clearvars; clear memoize; % clear memoize required for memoization
addpath(genpath('./'));

date.day            = 17;
date.mon            = 11;
date.year           = 2017;
timeStart           = '00:00:00.000';
timeEnd             = '05:00:00.000';
timeDataXRay        = TimeExt(timeStart, timeEnd, 100, date, 1e7); % change refreshSunMoonInfluenceTime to real number
iterationNumber     = 1;
secondInOneMinute   = 60;

%{'ukf', 'srukf', 'cdkf', 'srcdkf', 'ckf', 'sckf', 'fdckf', 'cqkf', 'ghqf', 'sghqf', 'ekf', 'pf', 'sppf', 'gspf', 'gmsppf'};
filterTypes = {'srukf'};

b_det   = 0.1; % Detector Background Rate. [photon*sec^-1]
b_diff  = 0.1; % Diffuse X-ray Background. [photon*cm^-2*sec^-1]
b_cosm  = 5; % Net Cosmic Ray Background. [photon*cm^-2*sec^-1]
x_ray_source_count      = 4;
time_bucket           = 1; % [sec]
detector_area         = 1; % [m^2]
background_photon_rate = b_det / (detector_area*1e4) + b_diff + b_cosm; % detector_area have to be converted to cm*cm

x_ray_sources   = load_x_ray_sources(x_ray_source_count);

initial_x_ray_orbit = load_initial_orbit();
initial_x_ray_orbit = initial_x_ray_orbit(1:6);

time_bucket_array = [1e1 1e2 1e3 1e4 1e5];
time_bucket_count = length(time_bucket_array);

dtoa_rms = zeros(x_ray_source_count, length(time_bucket_array));

for i = 1:time_bucket_count
    dtoa_cov = diag( x_ray_dtoa_covariance(x_ray_sources, detector_area, time_bucket_array(i), background_photon_rate) );
    dtoa_rms(:, i) = sqrt(dtoa_cov);
end

for i = 1:x_ray_source_count
    figure();
    loglog(time_bucket_array, dtoa_rms(i, :), 'LineWidth', 2);
    grid on;
    hold on;
    title(['1-sigma pulse DTOA measurement accurracy (sec) for ', x_ray_sources(i).Name]);
    xlabel('A*dt, m*sec');
    ylabel('1-sigma pulse DTOA, sec');
end

