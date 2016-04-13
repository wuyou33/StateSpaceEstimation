function [ condition ] = initialCondition( )
%%  Build initial condition for integrated INS and SNS navigation system

%%
    initialDistance   = [17; 10; 15];
    initialVelocity   = [0.01; 0.01; 0.01];
    initialQuaternion = [1; 0; 0; 0];
    initialAccelerationBias = [0; 0; 0];
    initialGyroBias = [0; 0; 0];
    initialAccelerationScale = [0; 0; 0];
    initialGyroScale = [0; 0; 0];
    
    condition = [initialDistance; initialVelocity; initialQuaternion; initialAccelerationBias; initialGyroBias; initialAccelerationScale; initialGyroScale];
end

