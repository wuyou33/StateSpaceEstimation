function [ result ] = fitSolarSystemGravityModel( sampleTime, sampleNumber )
    % fitSolarSystemGravityModel.
    %
    % This function is to fit gravity model of solar sytem.
    % To calculate gravity model function need to calcualte the motion of the planets based on their initial velocity,
    % as well as the gravitational influence of other planets, moons of earth, jupiter & saturn, as well as halleys comet.
    %
    %   INPUT
    %       initial             array of initial condition (position and velocity in in ICRF-j2000 coordinate system, and mass) for all planetary objects;
    %       sampleTime          sample time, [sec];
    %       sampleNumber        number of simulated samples.
    %
    %   OUTPUT
    %       result              SolarSystemGravityModel which contains array of GravityBodyInfo struct for all planetary objects from input arguments.
    %
    narginchk(2, 2);
    
    G = 6.67e-11; % gravitational constant [m^3kg^-1*s^-2].
    m_in_km = 1e3;
    
    initial = loadAstroBodyData();
    
    totalPlanets = length(initial);
    astrObjects  = zeros(totalPlanets, 6);
    massList = zeros(totalPlanets, 1);
    bodyOut  = zeros(totalPlanets, sampleNumber, 3);
    
    for i = 1 : totalPlanets;
        body = initial(i);
        astrObjects(i, 1) = body.X;
        astrObjects(i, 2) = body.Y;
        astrObjects(i, 3) = body.Z;
        
        astrObjects(i, 4) = body.VX;
        astrObjects(i, 5) = body.VY;
        astrObjects(i, 6) = body.VZ;
        
        massList(i) = body.Mass;
    end
    
    for time = 1 : sampleNumber
        for i = 1 : totalPlanets
            pos_i       = [astrObjects(i, 1), astrObjects(i, 2), astrObjects(i, 3)];
            vel_i       = [astrObjects(i, 4), astrObjects(i, 5), astrObjects(i, 6)];
            mass_i      = massList(i);
            
            gravList_x = zeros(totalPlanets);
            gravList_y = zeros(totalPlanets);
            gravList_z = zeros(totalPlanets);
            
            % next loop compares all of the planets / bodies (except 'i') with planet / body 'i'.
            % The stepping is wierd (a+1) as it is only necessary to calculate the grav vector from i to j.
            % Once calculated, it is easy to reverse the calculated vector (multiply by -1) to get the j to i vector.
            
            for j = (i+1) : totalPlanets
                diff_j_to_i = astrObjects(j, 1:3) - pos_i;
                pos_diff_mag = sqrt(diff_j_to_i(1)*diff_j_to_i(1) + diff_j_to_i(2)*diff_j_to_i(2) + diff_j_to_i(3)*diff_j_to_i(3));
                
                grav = diff_j_to_i / pos_diff_mag * ( (mass_i * massList(j) * G) / ((pos_diff_mag * m_in_km)^2) );
                
                % gravity between i and j
                gravList_x(i, j) = -grav(1);
                gravList_y(i, j) = -grav(2);
                gravList_z(i, j) = -grav(3);
                
                % gravity between j and i
                gravList_x(j, i) = grav(1);
                gravList_y(j, i) = grav(2);
                gravList_z(j, i) = grav(3);
            end
            
            % this turns the net grav value into an acceleration value [km/sec^2], by removing mass.
            totalGravAcceleration_i = [sum(gravList_x(:, i), 1), sum(gravList_y(:, i), 1), sum(gravList_z(:, i), 1)] / mass_i / m_in_km;
            
            astrObjects(i, 4 : 6) = vel_i + totalGravAcceleration_i * sampleTime;
            astrObjects(i, 1 : 3) = pos_i + (astrObjects(i, 4 : 6) * sampleTime);
        end
        
        for i = 1 : totalPlanets
            bodyOut(i, time, :) = astrObjects(i, 1:3);
        end
    end
    
    tmp = rvecrep(GravityBodyInfo(0, [0; 0; 0]), totalPlanets);
    for i = 1 : totalPlanets
        tmp(i) = GravityBodyInfo(massList(i), squeeze(bodyOut(i, :, :))');
    end
    
    result = SolarSystemGravityModel(tmp);
end
