function [S0] = s0_Nut( jd, nut) 
%Имя функции: s0_Nut 
%Функция рассчитывает истинное или среднее звездное время на 0ч UTC 
%Входные данные: 
%timeUTC - дата, на которую требуется рассчитать истинное или среднее звездное время 
%nut- признак (если nut= 0, то вычисляется среднее звездное время без учета нутации, иначе вычис-ляется 
%истинное звездное время) 
%Выходные данные: 
%S0 -  истинное или среднее звездное время на 0ч UTC  
jd2000 = 2451545; % 12h UTC 1 января 

d = jd - jd2000; 
 t = d / 36525.0;   % 36525 - юлианский период 100 лет  
    t2 = t * t; 
 h1 = 24110.54841;  
%h1=6.0*3600.0+41.0*60.0+50.54841;  
% h2 = 236.555367908 * d;     
   h2 = 8640184.812866 * t ;  
   h3 = 0.093104 * t2; 
   h4 = t2 * t * 6.2E-6; 
   if ( nut == 0)  
      na = 0; 
  else 
   na = utc_nut(t);  
  end;   
    s0_m = h1 + h2 + h3 - h4; 
    S0.s0_nut = s0_m + na; 
    S0.s0_m_mod = mod(s0_m, 86400); 
    s0_day = floor(s0_m / 86400); 
    S0.s0_m_hour = S0.s0_m_mod / 3600.0; 
    S0.s0_m_hour = floor(S0.s0_m_mod / 3600); 
    sec_min = S0.s0_m_mod - S0.s0_m_hour * 3600; 
    S0.s0_m_min = floor(sec_min / 60); 
    S0.s0_m_sec = sec_min - S0.s0_m_min * 60; 
    S0.s0_nut_mod = mod(S0.s0_nut, 86400); 
    s0_day = floor(S0.s0_nut / 86400); 
    S0.s0_nut_hour = S0.s0_nut_mod / 3600.0; 
    S0.s0_nut_hour = floor(S0.s0_nut_mod / 3600); 
    sec_min = S0.s0_nut_mod - S0.s0_nut_hour * 3600; 
    S0.s0_nut_min = floor(sec_min / 60); 
    S0.s0_nut_sec = sec_min - S0.s0_nut_min * 60; 
