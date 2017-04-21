%***PARA GUI***********
%17/06/2009
function show(dFFo_dat,co,s,M)

%*********línea para el umbral*********
lx=(1:M);
ffo_pos_pre=find(dFFo_dat<0);
ffo_pos=dFFo_dat;
ffo_pos(ffo_pos_pre)=0;

% co=abs(mean(ffo_pos));
% s=std(ffo_pos);

umbral=co;%*co/(s*0.5);
% umbral=6; %Determino el umbral a partir del ruido basal function show y spikes
ly(1:M)=umbral;
line(lx,ly,'Color','r')
%**********************************

%*********líneas de espigas*********
mx=max(dFFo_dat);
%Espigas

spikes_FFo=spikes(dFFo_dat,co,s);

plot(spikes_FFo*mx*1.5,'Color','g')