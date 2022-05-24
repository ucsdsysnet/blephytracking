function [amp,e,phi,I,Q,IQO,IQI,f0,phi_off,error,est_signal,fr,e_o] = BLE_Imperfection_Estimator_NAGD(x,seq,Fs,init_f0,init_e,init_phi,init_I,init_Q,init_amp,snr,n_partition)
% This function uses Nesterov Accelerated Gradient Descent to estimate
% hardware imperfection fingerprints of a Bluetooth signal
% Source: "Evaluating Physical-Layer BLE Location Tracking Attacks on Mobile Devices", IEEE SP 22

snr = 10^(snr/20);
err_thresh = max(0.4,1/(snr+1));

x2 = x(1:end-2*Fs/1e6);

t2 = (1:(length(x)))*(1/Fs);

seq = [1;0;seq];
% Creating the clean signal
est_signal_perfect2 = gfsk_modulate(seq,500e3,Fs).';
if n_partition ~= 1
    est_signal_perfect2 = (est_signal_perfect2(3.5*Fs/1e6+1:end-0.5*Fs/1e6));
else
    est_signal_perfect2 = (est_signal_perfect2(3.0*Fs/1e6+1:end-1*Fs/1e6));
end

fr = [];
I2=0;
Q2=0;
IQO2=0;
IQI2=0;
e2=0;
phi2=0;
phi_off2=0;
f02 = 0;
amp2 = 0;
e_o =[];
%freqsep2 = 0;
error2 = 0;
err = 0;

% Partitioning the samples to speed up the code and running over multiple
% random partitions to improve estimation accuracy and robustness
n = min(10,n_partition);
l = floor(length(x2)/n_partition);

for inter= 1:n
    r = randperm(length(est_signal_perfect2));
    r = r(1:l);
    est_signal_perfect = est_signal_perfect2(r);
    t = t2(r);
    x = x2(r);

    
    % Initialization
    e_new = init_e;
    phi_new = init_phi/180*pi;
    I_new = init_I;
    Q_new = init_Q;

    if inter == 1
        f0_new = init_f0;
    else
        f0_new = f0;
        init_f0 = f0;
    end
    w0_new = 2*pi*f0_new;
    phi_off_new = 36/360*2*pi;
    amp_new = init_amp;

    % Parameter update
    e = e_new;
    phi = phi_new;
    I = I_new;
    Q = Q_new;
    w0 = w0_new;
    phi_off = phi_off_new;
    amp = amp_new;
    %freqsep = freqsep_new;
    est_signal = ((amp-e)*(real(est_signal_perfect)*cos(phi)-imag(est_signal_perfect)*sin(phi))+...
        1i*(amp+e)*(imag(est_signal_perfect)*cos(phi)+real(est_signal_perfect)*sin(phi))+...
        I+1i*Q) .* exp(1i*(w0*t+phi_off));



    e_t = 0;
    phi_t = 0;
    I_t = 0;
    Q_t = 0;
    w0_t = 0;
    phi_off_t = 0;
    amp_t = 0;

    error_diff = 1;

    count = 0;
    round = 1;
    error = 1;
    
    while error_diff > 1e-7 && count<10e3%2e3
        count = count+1;
        
        % Setting the learning rate and momentum
        if error_diff > 1e-7
            lr = 1e-3;
            mom = 0.9;
        elseif error_diff < 1e-7
            lr = 1e-4;
            mom = 0.9;
        else
            lr = 1e-3;
            mom = 0.9;
        end


        % Momentum step
        %e = e_new - mom*e_t;
        phi = phi_new - mom*phi_t;
        I = I_new - mom*I_t;
        Q = Q_new - mom*Q_t;
        w0 = w0_new - mom*w0_t;
        phi_off = phi_off_new - mom*phi_off_t;
        %amp = amp_new - mom*amp_t;
        %freqsep = freqsep_new - mom*freqsep_t;

        % Re-initialization
        if count > round*2e2 && error(end) > err_thresh
            if floor(round/2)*2 == round-1
                round = round+1;
                f0 = init_f0-floor(round/2)*1.5e3;
                w0_new = 2*pi*f0;
            else
                round = round+1;
                f0 = init_f0+(round/2)*1.5e3;
                w0_new = 2*pi*f0;

            end
        end

        % Computing updated created imperfect signal
        Imag_part = ((amp-e)*(real(est_signal_perfect)*cos(phi)-imag(est_signal_perfect)*sin(phi))+I).*sin(w0*t+phi_off) +...
            ((amp+e)*(imag(est_signal_perfect)*cos(phi)+...
            real(est_signal_perfect)*sin(phi))+Q).*cos(w0*t+phi_off);

        Real_part = ((amp-e)*(real(est_signal_perfect)*cos(phi)-imag(est_signal_perfect)*sin(phi))+I).*cos(w0*t+phi_off) -...
            ((amp+e)*(imag(est_signal_perfect)*cos(phi)+...
            real(est_signal_perfect)*sin(phi))+Q).*sin(w0*t+phi_off);

        % Gradient descent calculation and update
        
        
        % freqsep_d = ((amp-e)*(*cos(phi)+imag(est_signal_perfect)*sin(phi))+I).*sin(w0*t+phi_off) +...
        %         ((amp+e)*(imag(est_signal_perfect)*cos(phi)+...
        %         real(est_signal_perfect)*sin(phi))+Q).*cos(w0*t+phi_off);

        e_d =  - mean((imag(x.')-Imag_part).*((-(real(est_signal_perfect)*cos(phi)-imag(est_signal_perfect)*sin(phi))).*sin(w0*t+phi_off)+...
            (imag(est_signal_perfect)*cos(phi) +...
            real(est_signal_perfect)*sin(phi)).*cos(w0*t+phi_off)));

        phi_d = -mean((imag(x.')-Imag_part).*((amp-e)*(-real(est_signal_perfect)*sin(phi)+imag(est_signal_perfect)*cos(phi)).*sin(w0*t+phi_off)+...
            ((amp+e)*(real(est_signal_perfect)*cos(phi) -...
            imag(est_signal_perfect)*sin(phi))).*cos(w0*t+phi_off)));

        I_d = -mean((imag(x.')-Imag_part).*(sin(w0*t+phi_off)));

        Q_d = -mean((imag(x.')-Imag_part).*(cos(w0*t+phi_off)));


        w0_d = -mean((imag(x.')-Imag_part).*(t.*((amp-e)*(real(est_signal_perfect)*cos(phi)-imag(est_signal_perfect)*sin(phi))+I).*cos(w0*t+phi_off) - ...
            t.*((amp+e)*(imag(est_signal_perfect)*cos(phi)+...
            real(est_signal_perfect)*sin(phi))+Q).*sin(w0*t+phi_off)));

        phi_off_d = -mean((imag(x.')-Imag_part).*(((amp-e)*(real(est_signal_perfect)*cos(phi)-imag(est_signal_perfect)*sin(phi))+I).*cos(w0*t+phi_off) - ...
            ((amp+e)*(imag(est_signal_perfect)*cos(phi)+...
            real(est_signal_perfect)*sin(phi))+Q).*sin(w0*t+phi_off)));

        amp_d = - mean((imag(x.')-Imag_part).* (((1)*(real(est_signal_perfect)*cos(phi)-imag(est_signal_perfect)*sin(phi))+I).*sin(w0*t+phi_off) +...
            ((1)*(imag(est_signal_perfect)*cos(phi)+...
            real(est_signal_perfect)*sin(phi))+Q).*cos(w0*t+phi_off)));


        e_d = e_d - mean((real(x.')-Real_part).*((-(real(est_signal_perfect)*cos(phi)-imag(est_signal_perfect)*sin(phi))).*cos(w0*t+phi_off)-...
            (imag(est_signal_perfect)*cos(phi) +...
            real(est_signal_perfect)*sin(phi)).*sin(w0*t+phi_off)));

        phi_d = phi_d - mean((real(x.')-Real_part).*((amp-e)*(-real(est_signal_perfect)*sin(phi)+imag(est_signal_perfect)*cos(phi)).*cos(w0*t+phi_off)-...
            ((amp+e)*(real(est_signal_perfect)*cos(phi) -...
            imag(est_signal_perfect)*sin(phi))).*sin(w0*t+phi_off)));

        I_d = I_d - mean((real(x.')-Real_part).*(cos(w0*t+phi_off)));

        Q_d = Q_d + mean((real(x.')-Real_part).*(sin(w0*t+phi_off)));


        w0_d = w0_d - mean((real(x.')-Real_part).*(-t.*((amp-e)*(real(est_signal_perfect)*cos(phi)-imag(est_signal_perfect)*sin(phi))+I).*sin(w0*t+phi_off) - ...
            t.*((amp+e)*(imag(est_signal_perfect)*cos(phi)+...
            real(est_signal_perfect)*sin(phi))+Q).*cos(w0*t+phi_off)));

        phi_off_d = phi_off_d - mean((real(x.')-Real_part).*(-((amp-e)*(real(est_signal_perfect)*cos(phi)-imag(est_signal_perfect)*sin(phi))+I).*sin(w0*t+phi_off) - ...
            ((amp+e)*(imag(est_signal_perfect)*cos(phi)+...
            real(est_signal_perfect)*sin(phi))+Q).*cos(w0*t+phi_off)));

        amp_d = amp_d - mean((real(x.')-Real_part).* (((1)*(real(est_signal_perfect)*cos(phi)-imag(est_signal_perfect)*sin(phi))+I).*cos(w0*t+phi_off) -...
            ((1)*(imag(est_signal_perfect)*cos(phi)+...
            real(est_signal_perfect)*sin(phi))+Q).*sin(w0*t+phi_off)));
        
        if error_diff >1e-7
            e_t = mom*e_t + lr*e_d;
            phi_t = mom*phi_t + lr*phi_d;
            I_t = mom*I_t + lr/1*I_d;
            Q_t = mom*Q_t + lr/1*Q_d;
            w0_t = mom*w0_t + 1e8*lr*w0_d;
            phi_off_t = mom*phi_off_t + 10*lr*phi_off_d;
            amp_t = mom*amp_t + lr*amp_d;
            %freqsep_t = mom*freqsep_t + lr*freqsep_d;
        else
            e_t = mom*e_t + lr*e_d;
            phi_t = mom*phi_t + lr*phi_d;
            I_t = mom*I_t + lr/1*I_d;
            Q_t = mom*Q_t + lr/1*Q_d;
            %w0_t = mom*w0_t + 1e8*lr*w0_d;
            phi_off_t = mom*phi_off_t + 10*lr*phi_off_d;
            amp_t = mom*amp_t + lr*amp_d;
        end

        %amp_new = amp + lr * mean((imag(x.')-Imag_part).* Imag_part / amp);



        e_new = e_new-e_t;
        phi_new = phi_new-phi_t;
        I_new = I_new-I_t;
        Q_new = Q_new-Q_t;
        w0_new = w0_new-w0_t;
        phi_off_new = phi_off_new-phi_off_t;
        amp_new = amp_new-amp_t;
        %freqsep_new = freqsep_new-freqsep_t;
        
        % updating parameters
        e = e_new;
        phi = phi_new;
        I = I_new;
        Q = Q_new;
        w0 = w0_new;
        phi_off = phi_off_new;
        amp = amp_new;
        %freqsep = freqsep_new;

        est_signal = ((amp-e)*(real(est_signal_perfect)*cos(phi)-imag(est_signal_perfect)*sin(phi))+...
        1i*(amp+e)*(imag(est_signal_perfect)*cos(phi)+real(est_signal_perfect)*sin(phi))+...
        I+1i*Q) .* exp(1i*(w0*t+phi_off));
        
        % computing the error
        %error = [error,mean(abs(est_signal.' - x)./abs(x))/2];
        error = [error,mean(norm(est_signal.'-x).^2)/mean(norm(x).^2)/2];

        if length(error) > 1 && error(end) < err_thresh
            error_diff = abs(error(end)-error(end-1));
        end
    end
    err = max([err,error(end)]);

    signal = x.*exp(-1j*(w0*t'+phi_off));

    try
        ell = fit_ellipse(real(signal),3*imag(signal));

        IQO = sqrt((ell.X0/ell.a)^2+(ell.Y0/ell.b)^2);
        IQI = ell.a/ell.b*3;
        flag = 1;
    catch

        warning('Ill ellipse');
        flag = 0;
    end

    if flag == 1
        f0 = w0/(2*pi);
        phi_off = phi_off*(360/(2*pi));
        phi = phi *(360/(2*pi));
        fr = [fr;f0];
        I2 = I2+I;
        Q2= Q2+Q;
        IQO2 = IQO2+IQO;
        IQI2 = IQI2+IQI;
        %e2 = e2+e;
        e2 = e2+e;
        phi_off2 = phi_off2+phi_off;
        phi2 = phi2+phi;
        f02 = f02+f0;
        amp2 = amp2+amp;
        error2 = error2+err;
    else
        disp('Ellipse was not found.');
    end
end
I = I2/n;
I = -I/amp;
Q = Q2/n;
Q = Q/amp;
IQO = IQO2/n;
IQI = IQI2/n;
phi=phi2/n;
e=e2/n;
phi_off=phi_off2/n;
f0 = f02/n;
amp = amp2/n;
error = err;
