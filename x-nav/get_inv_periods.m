function [ inv_periods ] = get_inv_periods( x_ray_sources )
    % get_inv_periods. Calculate following expression 2*pi / Period for every X-Ray Source.
    %
    %   [ inv_periods ] = get_inv_periods( x_ray_sources )
    %
    %   INPUT
    %       x_ray_sources   array of the X-Ray sources, every item should be instance of the <X_RaySource>.
    %
    %   OUTPUT
    %       inv_periods    result of 2*pi / Period.
    %
    dimension = length(x_ray_sources);
    inv_periods = zeros(1, dimension);
    
    for i = 1 : dimension
        x = x_ray_sources(i);
        inv_periods(i) = x.TwoPiOnPeriod;
    end
end
