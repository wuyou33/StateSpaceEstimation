classdef SatellitePhaseSpace < handle
    properties (Access = private)
        trajectory;
        velocity;
        rotation;
        iteration;
        capacity;
    end
    
    properties (Dependent)
        Trajectory;
        Velocity;
        Rotation;
        TrajectoryX;
        TrajectoryY;
        TrajectoryZ;
        VelocityX;
        VelocityY;
        VelocityZ;
        RotationO;
        RotationI;
        RotationJ;
        RotationK;
    end
    
    methods 
        function obj = SatellitePhaseSpace(phaseStateArray, capacity)
            obj.trajectory = zeros(3, capacity);
            obj.velocity = zeros(3,capacity);
            obj.rotation = zeros(4,capacity);
            
            [~, iterationNumber] = size(phaseStateArray);
            
            obj.trajectory(:, 1: iterationNumber) = phaseStateArray(1:3,:);
            obj.velocity(:, 1: iterationNumber) = phaseStateArray(4:6,:);
            obj.rotation(:, 1: iterationNumber) = phaseStateArray(7:10,:);
            
            obj.capacity = capacity;
            obj.iteration = iterationNumber;
        end
        
        function val = get.TrajectoryX(obj)
            val = obj.trajectory(1,:);
        end
        
        function val = get.TrajectoryY(obj)
            val = obj.trajectory(2,:);
        end
        
        function val = get.TrajectoryZ(obj)
            val = obj.trajectory(3,:);
        end
        
        function val = get.Trajectory(obj)
            val = obj.trajectory;
        end
        
        function val = get.Velocity(obj)
            val = obj.velocity;
        end
        
        function val = get.VelocityX(obj)
            val = obj.velocity(1,:);
        end
        
        function val = get.VelocityY(obj)
            val = obj.velocity(2,:);
        end
        
        function val = get.VelocityZ(obj)
            val = obj.velocity(3,:);
        end
        
        function val = get.Rotation(obj)
            val = obj.rotation;
        end
        
        function val = get.RotationO(obj)
            val = obj.rotation(1,:);
        end
        
        function val = get.RotationI(obj)
            val = obj.rotation(2,:);
        end
        
        function val = get.RotationJ(obj)
            val = obj.rotation(3,:);
        end
        
        function val = get.RotationK(obj)
            val = obj.rotation(4,:);
        end
        
        function val = GetPhaseState(obj, index)
            val = [obj.trajectory(:, index); obj.velocity(:, index); obj.rotation(:, index)];
        end
                
        function obj = AddPhaseState(obj, state, index)
            obj.trajectory(:, index) = state(1:3);
            obj.velocity(:, index) = state(4:6);
            obj.rotation(:, index) = state(7:10);            
        end
    end
    
end

