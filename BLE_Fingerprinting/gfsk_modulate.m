function y  = gfsk_modulate(x,freqsep,Fs)
% This function generates a BLE signal.
% Fs: sampling rate (Hz).
% freqsep: separation frequency of 0 and 1 bits.

nsample = Fs/1e6;
t = (1:(nsample*length(x)))*(1/Fs);

gamma_fsk = zeros(1,length(t));
for i=1:length(x)
    gamma_fsk((((i-1)*nsample)+1):(i*nsample)) = ((x(i)*2)-1);
end

gaussFilter = gaussdesign(0.3, 3, nsample); 
gamma_gfsk = filter(gaussFilter, 1, gamma_fsk);
gfsk_phase = (freqsep/Fs)*pi*cumtrapz(gamma_gfsk);
y = exp(1i*gfsk_phase);
y = y.';

end