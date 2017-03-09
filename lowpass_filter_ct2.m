function x = lowpass_filter_ct2(x)
% Use the filtfilt function so that we have zero phase lag.  Remember,
% that x is a matrix, but the filtfilt command can handle this.


% % Setup the options.
% if size(x,2)<1000
%     nlfilt=20;
% end
% if 1001<size(x,2)<2000
%     nlfilt=40;
% end
% if 2001<size(x,2)<3000
% nlfilt = 60; 
% end
%LC2014 Funciona bien para Fo de hasta 1000 frames
% options.nlfilt.value;		% filter length
% lpass = 0.25;
% % options.lpass.value;		% high pass in Hz.
% Fs = 4;
% 1/options.timeRes.value;		% 
nlfilt=40; %filtro de orden 40 la frecuencia de corte es casi cuadrada

[nrecordings len] = size(x);


%nlfilt = 5000;			% filter length.
if (nlfilt > len/3-1) % filt length must less than 1/3 of the data
    nlfilt = floor(len/3-1);
    if mod(nlfilt,2);
	nlfilt = nlfilt-1;
    end
end
%frecuencia de muestreo en milisegundos LCen14
normfreq_lpass = 0.15;%lpass/(Fs/2);	% 1 corresponds to Nyquist rate. La mitad de la 1/freq funciona
hfilt = fir1(nlfilt, normfreq_lpass, 'low');
xfilt = filtfilt(hfilt, 1, x')'; %zero phase digital filter LCen14
x = xfilt;
