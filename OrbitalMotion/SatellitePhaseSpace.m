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
        State;
        FullState;
    end
    
    methods 
        function obj = SatellitePhaseSpace(phaseStateArray, capacity)
            obj.trajectory = zeros(3, capacity);
            obj.velocity = zeros(3, capacity);
            obj.rotation = zeros(4, capacity);
            
            [~, iterationNumber] = size(phaseStateArray);
            
            obj.trajectory(:, 1:iterationNumber) = phaseStateArray(1:3, :);
            obj.velocity(:, 1:iterationNumber) = phaseStateArray(4:6, :);
            obj.rotation(:, 1:iterationNumber) = phaseStateArray(7:10, :);
            
            obj.capacity = capacity;
            obj.iteration = iterationNumber;
        end
        
        function val = get.State(this)
            val = [this.trajectory(1:3, :); this.velocity(1:3, :)];
        end
        
        function val = get.FullState(this)
            val = [this.trajectory(1:end, :); this.velocity(1:end, :); this.rotation(1:end, :)];
        end
        
        function val = get.TrajectoryX(this)
            val = this.trajectory(1, :);
        end
        
        function val = get.TrajectoryY(this)
            val = this.trajectory(2, :);
        end
        
        function val = get.TrajectoryZ(this)
            val = this.trajectory(3, :);
        end
        
        function val = get.Trajectory(this)
            val = this.trajectory;
        end
        
        function val = get.Velocity(this)
            val = this.velocity;
        end
        
        function val = get.VelocityX(this)
            val = this.velocity(1, :);
        end
        
        function val = get.VelocityY(this)
            val = this.velocity(2, :);
        end
        
        function val = get.VelocityZ(this)
            val = this.velocity(3, :);
        end
        
        function val = get.Rotation(this)
            val = this.rotation;
        end
        
        function val = get.RotationO(this)
            val = this.rotation(1, :);
        end
        
        function val = get.RotationI(this)
            val = this.rotation(2, :);
        end
        
        function val = get.RotationJ(this)
            val = this.rotation(3, :);
        end
        
        function val = get.RotationK(this)
            val = this.rotation(4, :);
        end
        
        function val = GetPhaseState(this, index)
            val = [this.trajectory(:, index); this.velocity(:, index); this.rotation(:, index)];
        end
                
        function AddPhaseState(this, state, index)
            this.trajectory(:, index) = state(1:3);
            this.velocity(:, index) = state(4:6);
            this.rotation(:, index) = state(7:10);            
        end
    end
    
end

