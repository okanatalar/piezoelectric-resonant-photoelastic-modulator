%% Processing for 2Vpp

spectrum_matrix = fft(channel_matrix_1_2Vpp,2e4*20,2);

[voltages_2Vpp,locations_2Vpp] = max(abs(spectrum_matrix),[],2);

phases_2Vpp = zeros(length(freq_points),1);
for k=1:size(spectrum_matrix,1)
    phases_2Vpp(k) = angle(spectrum_matrix(k,locations_2Vpp(k)));
end

figure;
plot(freq_points,voltages_2Vpp)

figure;
plot(freq_points,phases_2Vpp)

% Processing for 2Vpp Switch Open

spectrum_matrix = fft(channel_matrix_1_2Vpp_calib,2e4*20,2);%Oversampling rate of 20 for fft to not miss the peak

[voltages_2Vpp_calib,locations_2Vpp_calib] = max(abs(spectrum_matrix),[],2);

phases_2Vpp_calib = zeros(length(freq_points),1);
for k=1:size(spectrum_matrix,1)
    phases_2Vpp_calib(k) = angle(spectrum_matrix(k,locations_2Vpp_calib(k)));
end

figure;
plot(freq_points,voltages_2Vpp_calib)

figure;
plot(freq_points,phases_2Vpp_calib)


% s_11 Measurement at 2Vpp Using the Switch Open Measurements
voltage_2Vpp_calibrated = voltages_2Vpp./voltages_2Vpp_calib;
phase_2Vpp_calibrated = phases_2Vpp - phases_2Vpp_calib;

for e=1:length(phase_2Vpp_calibrated)
    if phase_2Vpp_calibrated(e) > 4
        phase_2Vpp_calibrated(e) = phase_2Vpp_calibrated(e) - 2*pi;
    elseif phase_2Vpp_calibrated(e) < -4
        phase_2Vpp_calibrated(e) = phase_2Vpp_calibrated(e) + 2*pi;
    end
end
        
figure;
plot(freq_points,voltage_2Vpp_calibrated)

figure;
plot(freq_points,phase_2Vpp_calibrated)

%Converting Voltage to Impedance
Z_2Vpp = 50*voltage_2Vpp_calibrated.*exp(1j*phase_2Vpp_calibrated)./(1 - voltage_2Vpp_calibrated.*exp(1j*phase_2Vpp_calibrated));

figure;
plot(freq_points,abs(Z_2Vpp))

%Converting Impedane to s_11
s_11_2Vpp = (Z_2Vpp - 50)./(Z_2Vpp + 50);

figure;
plot(freq_points,abs(s_11_2Vpp))


%% Processing for 4Vpp

spectrum_matrix = fft(channel_matrix_1_4Vpp,2e4*20,2);

[voltages_4Vpp,locations_4Vpp] = max(abs(spectrum_matrix),[],2);

phases_4Vpp = zeros(length(freq_points),1);
for k=1:size(spectrum_matrix,1)
    phases_4Vpp(k) = angle(spectrum_matrix(k,locations_4Vpp(k)));
end

figure;
plot(freq_points,voltages_4Vpp)

figure;
plot(freq_points,phases_4Vpp)


% s_11 Measurement at 4Vpp Using the Switch Open Measurements
voltage_4Vpp_calibrated = voltages_4Vpp./(2*voltages_2Vpp_calib);
phase_4Vpp_calibrated = phases_4Vpp - phases_2Vpp_calib;

for e=1:length(phase_4Vpp_calibrated)
    if phase_4Vpp_calibrated(e) > 4
        phase_4Vpp_calibrated(e) = phase_4Vpp_calibrated(e) - 2*pi;
    elseif phase_4Vpp_calibrated(e) < -4
        phase_4Vpp_calibrated(e) = phase_4Vpp_calibrated(e) + 2*pi;
    end
end
        
figure;
plot(freq_points,voltage_4Vpp_calibrated)

figure;
plot(freq_points,phase_4Vpp_calibrated)

%Converting Voltage to Impedance
Z_4Vpp = 50*voltage_4Vpp_calibrated.*exp(1j*phase_4Vpp_calibrated)./(1 - voltage_4Vpp_calibrated.*exp(1j*phase_4Vpp_calibrated));

figure;
plot(freq_points,abs(Z_4Vpp))

%Converting Impedane to s_11
s_11_4Vpp = (Z_4Vpp - 50)./(Z_4Vpp + 50);

figure;
plot(freq_points,abs(s_11_4Vpp))

%% Processing for 8Vpp

spectrum_matrix = fft(channel_matrix_1_8Vpp,2e4*20,2);

[voltages_8Vpp,locations_8Vpp] = max(abs(spectrum_matrix),[],2);

phases_8Vpp = zeros(length(freq_points),1);
for k=1:size(spectrum_matrix,1)
    phases_8Vpp(k) = angle(spectrum_matrix(k,locations_8Vpp(k)));
end

figure;
plot(freq_points,voltages_8Vpp)

figure;
plot(freq_points,phases_8Vpp)


% s_11 Measurement at 4Vpp Using the Switch Open Measurements
voltage_8Vpp_calibrated = voltages_8Vpp./(4*voltages_2Vpp_calib);
phase_8Vpp_calibrated = phases_8Vpp - phases_2Vpp_calib;

for e=1:length(phase_8Vpp_calibrated)
    if phase_8Vpp_calibrated(e) > 4
        phase_8Vpp_calibrated(e) = phase_8Vpp_calibrated(e) - 2*pi;
    elseif phase_8Vpp_calibrated(e) < -4
        phase_8Vpp_calibrated(e) = phase_8Vpp_calibrated(e) + 2*pi;
    end
end
        
figure;
plot(freq_points,voltage_8Vpp_calibrated)

figure;
plot(freq_points,phase_8Vpp_calibrated)

%Converting Voltage to Impedance
Z_8Vpp = 50*voltage_8Vpp_calibrated.*exp(1j*phase_8Vpp_calibrated)./(1 - voltage_8Vpp_calibrated.*exp(1j*phase_8Vpp_calibrated));

figure;
plot(freq_points,abs(Z_8Vpp))

%Converting Impedane to s_11
s_11_8Vpp = (Z_8Vpp - 50)./(Z_8Vpp + 50);

figure;
plot(freq_points,abs(s_11_8Vpp))

%% Processing for 12Vpp

spectrum_matrix = fft(channel_matrix_1_12Vpp,2e4*20,2);

[voltages_12Vpp,locations_12Vpp] = max(abs(spectrum_matrix),[],2);

phases_12Vpp = zeros(length(freq_points),1);
for k=1:size(spectrum_matrix,1)
    phases_12Vpp(k) = angle(spectrum_matrix(k,locations_12Vpp(k)));
end

figure;
plot(freq_points,voltages_12Vpp)

figure;
plot(freq_points,phases_12Vpp)


% s_11 Measurement at 4Vpp Using the Switch Open Measurements
voltage_12Vpp_calibrated = voltages_12Vpp./(6*voltages_2Vpp_calib);
phase_12Vpp_calibrated = phases_12Vpp - phases_2Vpp_calib;

for e=1:length(phase_12Vpp_calibrated)
    if phase_12Vpp_calibrated(e) > 4
        phase_12Vpp_calibrated(e) = phase_12Vpp_calibrated(e) - 2*pi;
    elseif phase_12Vpp_calibrated(e) < -4
        phase_12Vpp_calibrated(e) = phase_12Vpp_calibrated(e) + 2*pi;
    end
end
        
figure;
plot(freq_points,voltage_12Vpp_calibrated)

figure;
plot(freq_points,phase_12Vpp_calibrated)

%Converting Voltage to Impedance
Z_12Vpp = 50*voltage_12Vpp_calibrated.*exp(1j*phase_12Vpp_calibrated)./(1 - voltage_12Vpp_calibrated.*exp(1j*phase_12Vpp_calibrated));

figure;
plot(freq_points,abs(Z_12Vpp))

%Converting Impedane to s_11
s_11_12Vpp = (Z_12Vpp - 50)./(Z_12Vpp + 50);

figure;
plot(freq_points,abs(s_11_12Vpp))