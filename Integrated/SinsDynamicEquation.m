function [ dState ] = SinsDynamicEquation( ~, ... time
    state, ...
    accelerometerBiasDriftRate, ...    
    gyroBiasDriftRate, ...    
    acceleration, ...
    angularVelocity, ...
    quaternion)
%%
%   state(1:3) - position error
%   state(4:6) - velocity error
%   state(7:10) - quaternion error
%   state(11:13) - acceleration bias error
%   state(14:16) - gyro bias error
%   state(17:19) - acceleration scale factor error
%   state(20:22) - gyro scale factor error

%   acceleration and gyro bias error describing following model 
%       Brownian motion (sometimes called arithmetic Brownian motion or 
%       generalized Wiener process)):
%       dx/dt = u(t)dt+V(t)dW;
%%       
    dState = zeros(1, 22);
    
    q0 = state(7);
    q1 = state(8);
    q2 = state(9);
    q3 = state(10);
    
    dState(1:3) = state(4:6);
    
    dState(4:6) = ( eye(3) - quat2rotm(state(7:10)) ) * acceleration + ...
        quatrotate( quaternion, (state(17:19)'.* acceleration)' )' + ...
        quatrotate( quaternion, state(11:13) )';
           
    dState(7:10) = -0.5* [-q1 -q2 -q3; q0 -q3 q2; q3 q0 -q1; -q2 q1 q0] * ( ...
        quatrotate(quaternion, state(14:16)) + ...
        quatrotate(quaternion, (state(20:22)'.*angularVelocity)'))';
    
    dState(11:13) = accelerometerBiasDriftRate'.*state(11:13);
    dState(14:16) = gyroBiasDriftRate'.*state(14:16);
    dState(17:19) = zeros(1,3);
    dState(20:22) = zeros(1,3);
end

