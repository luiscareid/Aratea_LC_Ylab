%Obtiene las espigas a partir de un valor umbral en la funcion de
%entrada dFFo_dat:funcion de entrada/ co:numero de veces la desviacion
%estandar para el umbral/ s:desviacion estandar de la funcion de entrada

function spikes_FFo=spikes(dFFo_dat,co,s)

ffo_pos_pre=find(dFFo_dat<0);
ffo_pos=dFFo_dat;
ffo_pos(ffo_pos_pre)=0;

co=abs(mean(ffo_pos));
s=std(ffo_pos);

umbral=100;%*co/(s*0.5);
% umbral=5; %Determino el umbral a partir del ruido basal function show y spikes

spikes_FFo=dFFo_dat>=umbral;
spikes_FFo=spikes_FFo*1;
