function [ objects ] = load_astro_body_data()
    % loadAstroBodyData. Load initial conditions of all bodies of Solar system which are considered during simulation.
    %
    %   OUTPUT:
    %       objects         array which contains instances of <AstroBody> objects.
    %
    
    % load data from files which contains corrdinates (in ICRF-j2000 coordinate system), velocity (in ICRF-j2000 coordinate system), mass.
    load('solar_system.mat');
    
    count = size(planets_au, 1);
    objects = row_vector_replicate(AstroBody(zeros(3, 1), zeros(3, 1), 0), count);
    
    for i = 1 : count
        bodyData = planets_au(i,:);
        objects(i) = AstroBody(bodyData(1:3), bodyData(4:6), bodyData(7));
    end
end
