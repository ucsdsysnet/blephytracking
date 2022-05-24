function [fingerprint,bits] = BLE_Fingerprint(signal,snr,fs,preamble_detect,interp_fac,n_partition)
% The code takes a detected BLE signal and computes its fingerprints
% including CFO, I/Q offset and I/Q imbalance with several different methods

% center frequency of receiver
fc = 2.48e9;
% BLE channel
ch = 2.48e9;

% upsampling factor for the BLE signal
fs = fs*interp_fac;

% final output fingerprint
fingerprint = [];
    
% Normalizing the signal
signal = signal/mean(abs(signal));

% interpolating the signal
signal = interp(signal,interp_fac);

% Removing the channel offset and moving the signal fft to center
if interp_fac ~= 1
    signal_fft = fftshift(fft(signal));
    ls = length(signal);

    signal_fft_centered = zeros(ls,1);
    lcenter = floor(ls/2);
    lchannel = floor((ch-fc)/fs*ls+ls/2);
    lbandwidth = floor(2/(fs/1e6)*(ls-1)/2);
    signal_fft_centered(lcenter-lbandwidth:lcenter+lbandwidth) = signal_fft(lchannel-lbandwidth:lchannel+lbandwidth);
    signal = ifft(ifftshift(signal_fft_centered));
end
     
% Decoding the BLE signal        
[signal, signal_freq, bits] = BLE_Decoder(signal,fs,preamble_detect);

%Estimating CFO for initialization using preamble averaging method
pream = signal_freq(1:fs/1e6*8-1);
est_cfo = mean(pream(1:end));
est_cfo2 = est_cfo;
if abs(est_cfo2) > 100e3
    est_cfo2 = 0;
end

%Running the fingerprinting code
[amp,epsilon,phi,I,Q,IQO,IQI,f0,phi_off,error,~,~] = BLE_Imperfection_Estimator_NAGD(signal,bits,fs,est_cfo2,0,0,0,0,1,snr,n_partition);
 
tt = 0:1/fs:length(signal)/fs-1/fs;
signal = signal.*exp(-1j*(2*pi*f0*tt'+phi_off/(360/(2*pi))));
            

try
    sig = signal(randperm(length(signal)));
    ell = fit_ellipse(-real(sig)/amp*5,3*imag(sig)/amp);
    flag = 1;
catch
    warning('Ill ellipse');
    flag = 0;
end
            
               
            
angsig = angle(signal);
spl = 8;
quar = zeros(1,spl);
for sp = 1:spl/2
    e1 = (((angsig>((sp-1)*2*pi/spl))+(angsig<(sp*2*pi/spl)))==2);
    quar(1,sp) = mean(signal(e1));
    e1 = (((angsig>(-pi+(sp-1)*2*pi/spl))+(angsig<(-pi+sp*2*pi/spl)))==2);
    quar(1,sp+spl/2) = mean(signal(e1));
end


fingerprint_vec = [error,...
                    amp,...
                    f0,...
                    est_cfo,...
                    IQO,...
                    I,...
                    Q,...
                    sqrt(I^2+Q^2),...
                    IQI,...
                    epsilon,...
                    phi,...
                    ell.X0/ell.a,...
                    ell.Y0/ell.b,...
                    ell.X0_in/ell.a,...
                    ell.Y0_in/ell.b,...
                    sqrt((ell.X0/ell.a)^2+(ell.Y0/ell.b)^2),...
                    sqrt((ell.X0_in/ell.a)^2+(ell.Y0_in/ell.b)^2),...
                    ell.a*3/ell.b/5,...
                    ell.phi,...
                    real(mean(quar)),...
                    imag(mean(quar)),...
                    abs(mean(quar)),...
                    mean(real(signal)),...
                    mean(imag(signal)),...
                    abs(mean(real(signal))+1i*mean(imag(signal)))];

if length(fingerprint_vec) ~= 25
    flag=0;
end


if error(end) < 0.45 && flag == 1
    fingerprint = [fingerprint;fingerprint_vec];
end
