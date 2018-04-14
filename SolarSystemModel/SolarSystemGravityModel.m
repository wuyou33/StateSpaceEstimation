classdef SolarSystemGravityModel
    % <SolarSystemGravityModel>.
    % Allow to calculate Gravitational interaction in the Solar system.
    % contains array of <GravityBodyInfo> for all astro objects which sould be considered in gravitational interaction.
    
    properties (SetAccess = private)
        gravity_body_data;
    end
    
    methods
        function obj = SolarSystemGravityModel(data)
            obj.gravity_body_data = data;
        end
        
        function acc = eval_gravity_acceleration(this, current_time, obj_position_to_earth, obj_mass)
            % Calculate gravity acceleration ([km * sec^-2]) from all objects which are considered in the gravity model.
            if (nargin == 3)
                mass = 1; % The mass can be neglected
                warning('SolarSystemGravityModel :: eval_gravity_acceleration automatically assign mass to 1 kg');
            else
                mass = obj_mass;
            end
            
            G = 6.67e-11; % gravitational constant [m^3*kg^-1*s^-2].
            m_in_km = 1e3;
            
            total_astr_obj = length(this.gravity_body_data);
            grav_acc_list = zeros(3, total_astr_obj);
            
            % note: third element in GravityBodyData related to earth. need to rewrite.
            obj_position_to_ssb = obj_position_to_earth + this.gravity_body_data(3).Position(current_time);
            
            for i = 1 : total_astr_obj
                pos_diff =  this.gravity_body_data(i).Position(current_time) - obj_position_to_ssb;
                pos_diff_mag = sqrt(pos_diff(1)*pos_diff(1) + pos_diff(2)*pos_diff(2) + pos_diff(3)*pos_diff(3));
                unit_vect_diff = pos_diff / pos_diff_mag;
                
                grav_force_mag = ((mass * this.gravity_body_data(i).Mass * G) / ((pos_diff_mag * m_in_km)^2));
                
                grav_acc_list(:, i) = unit_vect_diff * grav_force_mag;
            end
            
            % the total gravitational force on body
            acc = sum(grav_acc_list, 2) / mass / m_in_km;
        end
        
    end
end
