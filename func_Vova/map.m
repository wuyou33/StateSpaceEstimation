function map(N) 
%Имя функции:map 
%Применена функция MatLab для внесения в графики орбитального движения изображения Земли 
load('topo.mat','topo','topomap1'); 
[x,y,z] = sphere(50); 
cla reset 
%axis square off 
props.AmbientStrength = 0.1; 
props.DiffuseStrength = 1; 
props.SpecularColorReflectance = .5;  
props.SpecularExponent = 20; 
props.SpecularStrength = 1; 
props.FaceColor= 'texture'; 
props.EdgeColor = 'none'; 
props.FaceLighting = 'phong'; 
props.Cdata = topo; 
surface(x*N,y*N,z*N,props); 
light('position',[-1 0 1]); 
light('position',[-1.5 0.5 -0.5], 'color', [.6 .2 .2]); 
view(3) ; 
grid on 
hold on