classdef SolarSystemGravityModel
    % SolarSystemGravityModel
    % Allow to calculate Gravitational interaction in the solar system.
    % contains array of <GravityBodyInfo> for all astro objects which sould be considered in gravitational interaction.
    
    properties (SetAccess = private)
        GravityBodyData;
    end
    
    methods
        function obj = SolarSystemGravityModel(data)
            narginchk(1, 1);
            
            obj.GravityBodyData = data;
        end
        
        function acc = EvalGravityAcceleration(this, sample, objPositionToEarth, objMass)
            % Calculate gravity acceleration ([km * sec^-2]) from all objects which are considered in the gravity model.
            if (nargin == 3)
                mass = 1; % The mass can be neglected
            else
                mass = objMass;
            end
            
            G = 6.67e-11; % gravitational constant [m^3kg^-1*s^-2].
            m_in_km = 1e3;
            
            totalAstrObj = length(this.GravityBodyData);
            gravAccList = zeros(3, totalAstrObj);
            
            % note: third element in GravityBodyData related to earth. need to rewrite.
            objPositionToSSB = objPositionToEarth + this.GravityBodyData(3).Position(:, sample);
            
            for i = 1 : totalAstrObj
                posDiff =  this.GravityBodyData(i).Position(:, sample) - objPositionToSSB;
                posDiffMag = sqrt(posDiff(1)*posDiff(1) + posDiff(2)*posDiff(2) + posDiff(3)*posDiff(3));
                unitVectDiff = posDiff / posDiffMag;
                
                gravForceMag = ((mass * this.GravityBodyData(i).Mass * G) / ((posDiffMag * m_in_km)^2));
                
                gravAccList(:, i) = unitVectDiff * gravForceMag;
            end
            
            % the total gravitational force on body
            acc = sum(gravAccList, 2) / mass / m_in_km;
        end
        
    end
end
