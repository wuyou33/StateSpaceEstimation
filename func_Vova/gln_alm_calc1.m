function [Pos_xyz,Vel_xyz] = gln_alm_calc1(ns,  n0, ti_current, alm_gln)      
%time_s0
% %Имя функции:gln_a_1 
% %Функция рассчитывает координаты и скорости спутников ГЛОНАСС в соответствии с интерфейс-ным  
% %контрольным документом ГЛОНАСС 
% %Входные данные: 
% %ns- номер спутника ,  
% %n0 - номер текущих суток внутри 4-х летнего периода (от ближайшего високосного года),  
% %ti_current -  текущее время обсервации от начала дня,  
% %time_s0 -  истинное звездное время в текущий момент обсервации, 
% %alm_gln - альманах спутников ГЛОНАСС 
% %Выходные данные: 
% %Структура satpos_xyz_gln  
% %satpos_xyz_gln.x - координата по оси x; 
% %satpos_xyz_gln.y- координата по оси y; 
% %satpos_xyz_gln.z- координата по оси z; 
% %satpos_xyz_gln.vx- скорость по оси x; 
% %satpos_xyz_gln.vy - скорость по оси y; 
% %satpos_xyz_gln.vz- скорость по оси z ; 
% % I_MID - Mean value of an inclination of a plane of orbit of a satellite 

global GL_MuE%[км^3/с^2] - константа гравитационного поля Земли
global GL_W_rot_Earth %[рад/с] - угловая скорость вращения Земли
global J2_0% - вторая зональная гармоника геопотенциала
global GL_Ae %[км] - экваториальный радиус Земли

  Tmid = 43200;%[с] - средний драконический период
   imid = 64.8;%[град] - среднее наклонение орбиты
   
   %3-шаг
   i =(imid/180 +  alm_gln.dI(ns))*pi;
   ecc_2 = alm_gln.E(ns)*alm_gln.E(ns);
   cos_i2 = cos(i)*cos(i);
   sin_i2 = sin(i)*sin(i);
   
   %1 - шаг`
   dtpr = (n0 - alm_gln.Na(ns)) * 86400 + (ti_current - alm_gln.tc(ns));
   
   %2 - шаг
   W = fix( dtpr./( Tmid + alm_gln.dT(ns) ) );
   
   %4 - шаг
   Tdr = Tmid + alm_gln.dT(ns) + (2.*W + 1) * alm_gln.dTT(ns);
   n = 2*pi./Tdr;
   
   %5 - шаг
   a_radius = semi_axis_1( Tdr, i, alm_gln.E(ns), alm_gln.omegan(ns)); 
%    alm_gln.a(ns) = a_radius;
   
   %6 - шаг
   p = a_radius .* (1 - ecc_2);
   %Долгота восхдящего узла
   Lamda =alm_gln.Lam(ns)*pi - ( GL_W_rot_Earth+3/2.*J2_0.*n.*( (GL_Ae./p).^2).*cos(i) ).*dtpr;
   %аргумент перигея
   Dae_p = GL_Ae./p;
   Koef = 3/4*J2_0.*n.*Dae_p.*Dae_p;
   w = alm_gln.omegan(ns)*pi - Koef.* (1-5*cos_i2).*dtpr;
   
   %7 
   E0 = -2* atan2(sqrt(1-alm_gln.E(ns)).*tan(w/2),sqrt(1+alm_gln.E(ns)));
   L1 = w + E0 - alm_gln.E(ns).* sin(E0);
   %8
   L_0  = L1 + n.*(dtpr - (Tmid + alm_gln.dT(ns)).*W - alm_gln.dTT(ns).*W.^2);
   
   %%9
   B = 3/2*J2_0*(GL_Ae./a_radius).^2;
   h = alm_gln.E(ns).*sin(w);
   l = alm_gln.E(ns).*cos(w);
   L(1,:) = L1;
   L(2,:) = L_0;
  
  D7_3 =  7.0 / 3.0;                  
  D7_4 = 7.0 / 4.0; 
  D7_6 = 7.0 / 6.0; 
  D7_24 = 7.0 / 24.0; 
  D49_72 = 49.0 / 72.0; 
  
  N_t=length(ti_current);
  
  deltaA=zeros(2,N_t);
  deltaH=zeros(2,N_t);
  deltal=zeros(2,N_t);
  deltaLam=zeros(2,N_t);
  deltaI=zeros(2,N_t);
  deltaL=zeros(2,N_t);
  
for j =1:2
  sin_l  = sin(L(j,:)); 
  sin_2l = sin(2.*L(j,:)); 
  sin_3l = sin(3.* L(j,:)); 
  sin_4l = sin(4.* L(j,:)); 
  cos_l  = cos(L(j,:)); 
  cos_2l = cos(2.* L(j,:)); 
  cos_3l = cos(3.* L(j,:)); 
  cos_4l = cos(4.* L(j,:));
  
   deltaA(j,:) = 2.*B.*(1-3/2*sin_i2).*(l.*cos_l + h.*sin_l)+ ...
               B.*sin_i2.*(1/2.*h.*sin_l-1/2.*l.*cos_l+  cos_2l+ ...
               7/2.*l.*cos_3l+ 7/2.*h.*sin_3l);
            
   deltaH(j,:) = B.*(1-3/2*sin_i2).*( sin_l+3/2.*l.*sin_2l-3/2.*h.*cos_2l)-...
               1/4.*B.*sin_i2.*( sin_l - D7_3.*sin_3l + 5.*l.*sin_2l -...
               17/2.*l.*sin_4l+17/2.*h.*cos_4l+h.*cos_2l)+...
               (-1/2.*B.*cos_i2.*l.*sin_2l);
            
   deltal(j,:) = B.*(1-3/2*sin_i2).*( cos_l+3/2.*l.*cos_2l+3/2.*h.*sin_2l)-...
               1/4.*B.*sin_i2.*( -cos_l - 7/3.*cos_3l - 5.*h.*sin_2l -...
               17/2.*l.*cos_4l-17/2.*h.*sin_4l+l.*cos_2l)+...
               (1/2.*B.*cos_i2.*h.*sin_2l);%???????
            
 deltaLam(j,:) = -B.*cos(i).*( 7/2.*l.*sin_l-5/2.*h.*cos_l-1/2.*sin_2l - ...
               7/6.*l.*sin_3l + 7/6.*h.*cos_3l);
   
   deltaI(j,:) = 1/2.*B.*cos(i).*sin(i).*( -l.*cos_l+h.*sin_l+cos_2l - ...
               7/3.*l.*cos_3l + 7/3.*h.*sin_3l);
                    %n*
   deltaL(j,:) = 2.0.*B.*(1-3/2*sin_i2).*...
               ( D7_4.*l.*sin_l - D7_4.*h.*cos_l) + 3.*B.*sin_i2.*... 
               ( - D7_24.*h.*cos_l - D7_24.*l.*sin_l -... 
               D49_72.*h.*cos_3l + D49_72.*l.*sin_3l +... 
               0.25.*sin_2l) + B.*cos_i2.*... 
               (3.5.*l.*sin_l - 2.5.*h.*cos_l -... 
               0.5.*sin_2l + D7_6.*l.*sin_3l + D7_6.*h.*cos_3l); %?????
end;
   LL = L_0;
   da = a_radius+deltaA(2,:).*a_radius-deltaA(1,:).*a_radius;  %!!!!!!!!!!!!!!!!!!!!!!!!!!
   dh = h + deltaH(2,:) - deltaH(1,:);
   dl = l + deltal(2,:) - deltal(1,:);
dLamda= Lamda + deltaLam(2,:) - deltaLam(1,:);
   di = i + deltaI(2,:) - deltaI(1,:);
   de = sqrt(dh.^2+dl.^2);
   dw = atan2(dh,dl);
   dL = LL + deltaL(2,:) - deltaL(1,:);  
   %kepler
   E_pp = dL - dw;
   E_npp = dL-dw+de.*sin(E_pp);
   
   for j=1:length(ti_current)
       while (abs( E_npp(j) - E_pp(j) )>=10^-9)
           E_pp(j) = E_npp(j);    
           E_npp(j) = dL(j)-dw(j)+de(j)*sin( E_pp(j) );
       end;
   end;
   %11 - истенная аномалия
   
   v = 2.*atan2( sqrt(1+de).*tan(E_npp/2),sqrt(1-de));   %истенная аномалия
   u = v+dw;                                           %аргумет широты
   
   %12 - истенные координаты
   p = da.*(1 - de.^2);                               %
   r = p./(1+de.*cos(v));
   
   v_r = sqrt(GL_MuE./p).*de.*sin(v);
   v_u = sqrt(GL_MuE./p).*(1+de.*cos(v));
   
   
   Pos_xyz(1,:) = r.*(cos(dLamda).*cos(u)- sin(dLamda).*sin(u).*cos(di));
   Pos_xyz(2,:) = r.*(sin(dLamda).*cos(u)+ cos(dLamda).*sin(u).*cos(di));
   Pos_xyz(3,:)= r.*sin(u).*sin(di);

   Vel_xyz(1,:) = v_r.*(cos(dLamda).*cos(u)-sin(dLamda).*sin(u).*cos(di))-...
                       v_u.*(cos(dLamda).*sin(u)+sin(dLamda).*cos(u).*cos(di))+...
                       GL_W_rot_Earth.*Pos_xyz(2,:);
   Vel_xyz(2,:) = v_r.*(sin(dLamda).*cos(u)+cos(dLamda).*sin(u).*cos(di))-...
                       v_u.*(sin(dLamda).*sin(u)-cos(dLamda).*cos(u).*cos(di))-...
                       GL_W_rot_Earth.*Pos_xyz(1,:); %???????
   Vel_xyz(3,:) = v_r.*sin(u).*sin(di) + v_u.*cos(u).*sin(di);