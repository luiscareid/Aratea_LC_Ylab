% FFo. Junio 2009
% Obtiene los transitorios de calcio de los puntos seleccionados LC
% El video *.tif tiene que estar en el mismo directorio
% 'video' es una cadena con el nombre del video, tiene que estar en el
% workspace
% El archivo de coordenadas 'coord' debe estar en el workspace, es un
% archivo nX2 donde n es el numero de celulas y las dos columnas
% corresponden a X,Y
%Inverti el orden en el algoritmo, primero abro la imagen y luego todos los
%centros para cada imagen, de esa forma cada imagen solo se abre unan vez.
%LC 04Dic09

function [FFo_dat]=FFo(coord,video)

h=4; %Mitad del cuadrado en pixeles

XY=coord;                            %selecciona centros
N=length(XY);                           %N centros
info=imfinfo(video);
M=length(info); %numero de frames
ajuste=info(1).Height; %altura de la imagen

%Fo/F
for i=1:M
    imagen=imread(video,i); %Abro solo una vez
    for j=1:N
    x1=XY(j,1);
    y1=XY(j,2); %Coordenadas obtenidas con Image J
   % y1=ajuste-y1; %Ajuste de coordenadas obtenidas con IDL    
 
    I=imcrop(imagen,[(x1-h) (y1-h) 2*h 2*h]);
    %promedio del area seleccionada.
    prom(i,j)=mean2(I);
    end
    
end

FFo_dat=prom;
