clc; close all; clearvars;

quaternion = [1;0;0;0];
N = 300;
w = randn(N, 3);

% quatFun = @(q, t, wb) -0.5 * [0 wb(1) wb(2) wb(3); -wb(1) 0 -wb(3) wb(2); -wb(2) wb(3) 0 -wb(1); -wb(3) -wb(2) wb(1) 0 ] * q;
quatFun = @(q, t, wb) -0.5 * quaternionMultiply(q, quaternionNormalize([0; wb]));

dt = 1e-1;
quat2 = quaternion;
quat1  = quaternion;
% quatmultiply
q1 = zeros(N, 4);
q2 = zeros(N, 4);
for i = 1: N
    wi = w(i, :);
    wnorm = norm(wi);
    co=cos(0.5*wnorm*dt);
    si=sin(0.5*wnorm*dt);
    
    n1 = wi(1)/wnorm;
    n2 = wi(2)/wnorm;
    n3 = wi(3)/wnorm;
    
    qw1 = n1*si;
    qw2 = n2*si;
    qw3 = n3*si;
    qw4 = co;
    
    Om=[ qw4  qw3 -qw2 qw1;
        -qw3  qw4  qw1 qw2;
         qw2 -qw1  qw4 qw3;
        -qw1 -qw2 -qw3 qw4];
    
    quat1 = Om * quat1;
    q1(i, :) = quat1;
    
    odeFun = @(t, q) quatFun(q, t, wi');
    [~, tmp] = ode45(odeFun, [(i-1)*dt; i*dt], quat2);
    quat2 = tmp(end, :)';
    q2(i, :) = quat2;
end

errors = q1 - q2;

figure();
plot(errors(:, 1)); grid on; hold on;

figure();
plot(errors(:, 2)); grid on; hold on;

figure();
plot(errors(:, 3)); grid on; hold on;

figure();
plot(errors(:, 4)); grid on; hold on;
