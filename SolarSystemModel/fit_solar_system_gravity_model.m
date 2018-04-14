function [ result ] = fit_solar_system_gravity_model( time_data )
    % fit_solar_system_gravity_model.
    %
    % This function is to fit gravity model of solar sytem.
    % To calculate gravity model function need to calcualte the motion of the planets based on their initial velocity,
    % as well as the gravitational influence of other planets, moons of earth, jupiter & saturn, as well as halleys comet.
    %
    %   [ result ] = fit_solar_system_gravity_model( time_data )
    %
    %   INPUT
    %       time_data   information about simulation time (start time, end time, simulation number, sample time, etc).
    %
    %   OUTPUT
    %       result  intance of the <SolarSystemGravityModel> which contains array of instances of the <GravityBodyInfo>
    %                   for all planetary objects from input arguments.
    %
    narginchk(1, 1);
    
    grav_const = 6.67e-11; % gravitational constant [m^3kg^-1*s^-2].
    m_in_km = 1e3;
    
    % assumptions: works correctly only when data in file solar_system.mat correspond to time from time_data.
    % at this moment need to sync these data manually. need to replace to call to ssd.jpl.nasa.gov/horizons.cgi (web interface which provides ephemeris).
    initial = load_astro_body_data();
    
    total_planets = length(initial);
    astr_objects  = zeros(total_planets, 6);
    mass_list = zeros(total_planets, 1);
    body_out  = zeros(total_planets, time_data.SimulationNumber, 3);
    
    for i = 1 : total_planets;
        body = initial(i);
        astr_objects(i, 1) = body.X;
        astr_objects(i, 2) = body.Y;
        astr_objects(i, 3) = body.Z;
        
        astr_objects(i, 4) = body.VX;
        astr_objects(i, 5) = body.VY;
        astr_objects(i, 6) = body.VZ;
        
        mass_list(i) = body.Mass;
    end
    
    for time = 1 : time_data.SimulationNumber
        for i = 1 : total_planets
            pos_i  = [astr_objects(i, 1), astr_objects(i, 2), astr_objects(i, 3)];
            vel_i  = [astr_objects(i, 4), astr_objects(i, 5), astr_objects(i, 6)];
            mass_i = mass_list(i);
            
            grav_list_x = zeros(total_planets);
            grav_list_y = zeros(total_planets);
            grav_list_z = zeros(total_planets);
            
            % next loop compares all of the planets / bodies (except 'i') with planet / body 'i'.
            % The stepping is wierd (a+1) as it is only necessary to calculate the grav vector from i to j.
            % Once calculated, it is easy to reverse the calculated vector (multiply by -1) to get the j to i vector.
            
            for j = (i+1) : total_planets
                diff_j_to_i = astr_objects(j, 1:3) - pos_i;
                pos_diff_mag = sqrt(diff_j_to_i(1)*diff_j_to_i(1) + diff_j_to_i(2)*diff_j_to_i(2) + diff_j_to_i(3)*diff_j_to_i(3));
                
                grav = diff_j_to_i / pos_diff_mag * ( (mass_i * mass_list(j) * grav_const) / ((pos_diff_mag * m_in_km)^2) );
                
                % gravity between i and j
                grav_list_x(i, j) = -grav(1);
                grav_list_y(i, j) = -grav(2);
                grav_list_z(i, j) = -grav(3);
                
                % gravity between j and i
                grav_list_x(j, i) = grav(1);
                grav_list_y(j, i) = grav(2);
                grav_list_z(j, i) = grav(3);
            end
            
            % this turns the net grav value into an acceleration value [km/sec^2], by removing mass.
            total_grav_acceleration_i = [sum(grav_list_x(:, i), 1), sum(grav_list_y(:, i), 1), sum(grav_list_z(:, i), 1)] / mass_i / m_in_km;
            
            astr_objects(i, 4 : 6) = vel_i + total_grav_acceleration_i * time_data.SampleTime;
            astr_objects(i, 1 : 3) = pos_i + (astr_objects(i, 4 : 6) * time_data.SampleTime);
        end
        
        for i = 1 : total_planets
            body_out(i, time, :) = astr_objects(i, 1:3);
        end
    end
    
    for i = 1 : total_planets
        tmp(i) = GravityBodyInfo(mass_list(i), squeeze(body_out(i, :, :))', time_data);
    end
    
    result = SolarSystemGravityModel(tmp);
end
