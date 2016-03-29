function jden = JD_epohi(epoha) 
%Имя: JD_epohi   
%Фукция JD_epohi(epoha) рассчитывает номер юлианского дня 
%до начала  года (epoha) на  12h, 0 день, январь. 
%Входные данные: epoha, размерность-год 
%Выходные данные: 
% jden- номер юлианского дня на  12h, 0 день, январь ( размерность -дни) 
  rk = mod(epoha,4); 
  if ( rk == 0 ) rk = 1.0; 
  else 
       rk = 2.0 - rk * 0.25; 
  end;  
   n100 = floor(epoha / 100); 
   n400 = floor(epoha / 400); 
   jden = (4712 + epoha) * 365.25 + n400 - n100 + rk; 
% fprintf('epoha=%d rk=%f jden=%6.2f \n', epoha, rk, jden); 
