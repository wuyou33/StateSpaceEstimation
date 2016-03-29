function y=S_fun(gmst,t) 
% расчёт угла поворота неинерциальной СК относительно инерциальной СК для заданного звёздного времени 
global GL_W_rot_Earth
y=gmst + GL_W_rot_Earth.*(t - 10800);