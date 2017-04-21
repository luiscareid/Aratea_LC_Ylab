% LC Enero 2010. A partir de una matriz de 1's y 0's quita todos los 1's
% consecutivos es decir, los que forman parte del burst. Con ello se
% facilita la seleccion de estados y la identificacion de los 
% ensambles neuronales en el LLE.

function [Spikes_trail]=trail(Spikes)

Sii=size(Spikes,1); %Numero de celulas
Sff=size(Spikes,2); %Numero de frames
n=0;                %Inicio de bandera
Spikes_trail=zeros(Sii,Sff);

for ii=1:Sii
    for jj=1:Sff
        if Spikes(ii,jj)==1
            n=n+1;
        end
        if Spikes(ii,jj)==0
            n=0;
        end        
        if n==1
            Spikes_trail(ii,jj)=1;
        end
    end;
    n=0; %Para evitar errores si el ultimo valor es 1
end;