function [ dState ] = SinsDynamicEquation( ~, state, acceleration, angularVelocity, quaternion)
%%
% SinsDynamicEquation - state transition function for integrated inertial navigation system (INS) and
% satellity navigation system (SNS). Describe ins dynamic error.
% 
%   state(1:3)      - position error
%   state(4:6)      - velocity error
%   state(7:10)     - quaternion error
%   state(11:13)    - acceleration bias error
%   state(14:16)    - gyro bias error
%   state(17:19)    - acceleration scale factor error
%   state(20:22)    - gyro scale factor error
% 
%   acceleration and gyro bias error describing following model 
%       Brownian motion (sometimes called arithmetic Brownian motion or 
%       generalized Wiener process)):
%       dx/dt = u(t)dt+V(t)dW;
% Equation (on the left side derivative by time):
%   dR = dV
%   dV = [I - C(Qpc)]*fp + C(Qbp)*[sf*fb + bf];
%   dQpc = -0.5*B*[C(Qbp)*(sw*wb + bw)];
%   where:
%       dR  - position error.
%       dV  - velocity error.
%       Qpc - quaternion roation from platform (p) frame to computed (c) frame.
%       fp  - acceleration resolved in platform (p) frame.
%       fb  - acceleration resolved in platform (b) frame.
%       Qbp - quaternion roation from body (b) frame to platform (p).
%       sf  - accelerometer scale factor.
%       bf  - accelerometer bias resolved in body (b) frame.
%       sw  - gyro scale factor.
%       wb  - angular velocity resolved in body (b) frame.
%       bw  - accelerometer bias resolved in body (b) frame.
%       B   - matrix from quaternion Qpc. (0.5*[-q1 -q2 -q3; q0 -q3 q2; q3 q0 -q1; -q2 q1 q0])
%   
%%  
    dState(1:3, 1) = state(4:6);
    
    % rotation matrix from body frame to platform frame
    body2Plat = quaternion2RotationMatrix(quaternion);
    
    % rotation matrix from platform frame to computer frame
    plat2Comp = quaternion2RotationMatrix(state(7:10));
    
    dState(4:6, 1) = (eye(3) - plat2Comp)*(body2Plat*acceleration) + body2Plat*(state(17:19).*acceleration + state(11:13));
    
    dState(7:10, 1) = -quaternion2BMatrix(state(7:10))*( body2Plat*(state(20:22).*angularVelocity + state(14:16)) );
    
    dState(11:13, 1) = [0; 0; 0];
    dState(14:16, 1) = [0; 0; 0];
    dState(17:19, 1) = [0; 0; 0];
    dState(20:22, 1) = [0; 0; 0];
end
