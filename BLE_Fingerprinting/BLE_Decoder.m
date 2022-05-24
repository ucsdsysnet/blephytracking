function [ble_signal, signal_freq, bits] = BLE_Decoder(signal,fs, preamble_detect)
% This function implements a very simple Bluetooth demodulator.
% More robust decoders can be used instead.
% The code assumes the input signal only consists of 1 Bluetooth signal.
% fs: sampling frequency (Hz)
% preamble_detect: If set to 1, the decoder also decodes the preamble. If
% set to zero, it assumes the Bluetooth signal starts from the first sample

%creating preamble
pream = [0,1,0,1,0,1,0,1,0,1,0];
preamble_signal = gfsk_modulate(pream,500e3,fs);
preamble_signal = preamble_signal(fs/1e6*2.5:end-fs/1e6*0.5);

signal_angle = unwrap(angle(preamble_signal));
slope = signal_angle(2:length(signal_angle))-signal_angle(1:length(signal_angle)-1);
preamble_freq = slope/(2*pi)*fs;
%figure; plot(preamble_freq);

% obtaining the instantaneous frequency of the signal
signal_angle = unwrap(angle(signal));
slope = signal_angle(3:length(signal_angle))-signal_angle(2:length(signal_angle)-1);
signal_freq = slope/(2*pi)*fs;
signal_freq = [signal_freq;0];
    
% Finding the preamble. 
% The code assumes the beginning of the packet has been detected almost
% accurately.
if preamble_detect == 0
    start_ind = 1;
else
    l = length(signal_freq);
    z = xcorr(signal_freq,preamble_freq);
    z = z(l+1:end);
    if length(z) > 20e-6*fs
        [~,start_ind] = max(abs(z(floor(2e-6*fs):floor(20e-6*fs))));
        %start_ind = start_ind+fs/1e6/2;
    else
        start_ind = 1;
    end
end

    signal = signal(start_ind:end);
    signal_freq = signal_freq(start_ind:end);
    

    signal_freq = signal_freq(1:floor(length(signal_freq)/(fs/1e6))*(fs/1e6));
    ble_signal = signal(1:floor(length(signal_freq)/(fs/1e6))*(fs/1e6));
    
    bits_freq = reshape(signal_freq,fs/1e6,length(signal_freq)/(fs/1e6))';
    bits = (mean(bits_freq(:,fs/1e6/2:fs/1e6/2+1),2)>0);
    %bits = (mean(bits_freq(:,40:60),2)>0);
    %bits = (mode(bits_freq(:,1:100)>0,2)>0);
    
end