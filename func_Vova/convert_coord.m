function [oxyz] = convert_coord(geogr_coord);
%функция преобразования географических координат в координаты ox0y0z0
%Входные данные:
%Структура geogr_coord
%geogr_coord.lat - широта  (градусы)
%geogr_coord.lon - долгота (градусы)
%geogr_coord.h   - высота
%
%Выходные данные:
%Структура oxyz
%oxyz.x 
%oxyz.y
%oxyz.z

eEarth = 0.01671123;        %екцентриситет земли
a = 6378245.0;              %большая полуось земли (метры)

e_2 = eEarth^2;
sin_lat2 = sin(geogr_coord.lat)*sin(geogr_coord.lat);
N = a/sqrt(1-e_2*sin_lat2); %радиус кривизны первого вертикала

oxyz.x = (N+geogr_coord.h)*cos(geogr_coord.lat)*cos(geogr_coord.lon);
oxyz.y = (N+geogr_coord.h)*cos(geogr_coord.lat)*sin(geogr_coord.lon);
oxyz.z = ((1-e_2)*N+geogr_coord.h)*sin(geogr_coord.lat);