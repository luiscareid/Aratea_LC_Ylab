function [ Spikes_fop ] = spikes_fast_oopsi( FFo,Freq,tau,smth )
%Dec 13 LC. Creado como una solucion alternativa al criterio de la derivada
%para encontrar las espigas de los transitorios de calcio
% set simulation metadata
V.dt    = 1/Freq;  % time step size
%initialize params
% P.a     = .2;    % observation scale
% P.b     = 0;    % observation bias
% tau     = 8;    % decay time constant
P.gam   = 1-V.dt/tau; % C(t) = gam*C(t-1)
% P.lam   = .1;  % firing rate = lam*dt
% P.sig   = 0.05;  % standard deviation of observation noise 
% F=FFo/max(FFo); %LC
F=FFo;
fr=smth;
window = ones(fr,1)/fr; 

 %C = filter(ones(1,fr)/fr,1,F); %filtra cada fr frames LC
 C=convn(F,window,'valid'); 

% fast oopsi
[Nhat Phat] = fast_oopsi(C,V,P);

Spikes_fop=Nhat;
end

