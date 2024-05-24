% Point 1 (Signal Generation)
t = 0:0.01:2*pi;
Amp = 10;
freq = 1;

% Generating a Rectangular Wave
rect_wave = square(2 * pi * freq * t);

% Generating a Complex Signal
amplitude_complex = 1.5;
complex_signal = amplitude_complex * exp(1j * 2 * pi * freq * t);

%---------------------------------------
% Point 2 (Addition - Sub - Multi)
Addition = rect_wave + real(complex_signal);
Subtraction = rect_wave - real(complex_signal);
Multi = rect_wave .* real(complex_signal);

figure;
subplot(4, 2, 1);
plot(t, rect_wave);
title('Rectangular Wave');
xlabel('Time');
ylabel('Amplitude');

subplot(4, 2, 2);
plot(t, real(complex_signal));
title('Real Part of Complex Signal');
xlabel('Time');
ylabel('Amplitude');

subplot(4, 2, 3);
plot(t, real(Addition));
title('Addition');
xlabel('Time');
ylabel('Amplitude');

subplot(4, 2, 4);
plot(t, real(Subtraction));
title('Subtraction');
xlabel('Time');
ylabel('Amplitude');

subplot(4, 2, 5);
plot(t, real(Multi));
title('Multiplication');
xlabel('Time');
ylabel('Amplitude');

%---------------------------------------
% Point 3 (Sampling of Signals)
frequency = 10; 
sampled_rect_wave = square(2 * pi * freq * (0:1/frequency:2));
sampled_complex_signal = real(complex_signal);

figure;
subplot(2, 2, 1);
t_sampled = 0:1/frequency:2;
plot(t_sampled, sampled_rect_wave, 'r', 'LineWidth', 1.5);
title('Sampled Rectangular Wave');
xlabel('Time (s)');
ylabel('Amplitude');

subplot(2, 2, 2);
stem(t_sampled, sampled_complex_signal(1:length(t_sampled)), 'b', 'LineWidth', 1.5);
title('Sampled Complex Signal');
xlabel('Time (s)');
ylabel('Amplitude');

%---------------------------------------
% Point 4 (Convolution & conv command comparison)
unit_step = ones(size(t));
unit_step(t < 0.5) = 0;

convolution_rect_wave = conv(rect_wave, unit_step);
convolution_complex_signal = conv(real(complex_signal), unit_step);

figure;
subplot(3, 3, 1);
plot(t, rect_wave);
title('Rectangular Wave');

subplot(3, 3, 3);
plot(t, real(complex_signal));
title('Real Part of Complex Signal');

subplot(3, 3, 5);
plot(convolution_rect_wave);
title('Convolution of Rectangular Wave');

subplot(3, 3, 2);
plot(convolution_complex_signal);
title('Convolution of Complex Signal');

subplot(3, 3, 4);
plot(t, rect_wave);
hold on;
plot(convolution_rect_wave(1:length(t)), 'r');
hold off;
title('Comparison - Rectangular Wave');

subplot(3, 3, 6);
plot(t, real(complex_signal));
hold on;
plot(convolution_complex_signal(1:length(t)), 'r');
hold off;
title('Comparison - Complex Signal');

% Point 5 (average filter)
noisy_signal_rect_wave = randn(1, length(t)); 
noisy_signal_complex_signal = randn(1, length(t));

filter_length = 10; 
filter_coefficients = ones(1, filter_length) / filter_length; 

filtered_rect_wave = filter(filter_coefficients, 1, noisy_signal_rect_wave);
filtered_complex_signal = filter(filter_coefficients, 1, noisy_signal_complex_signal);

figure;
subplot(2, 1, 1);
plot(noisy_signal_rect_wave);
title('Noisy Rectangular Wave');

subplot(2, 1, 2);
plot(filtered_rect_wave);
title('Filtered Rectangular Wave');

figure;
subplot(2, 1, 1);
plot(noisy_signal_complex_signal);
title('Noisy Complex Signal');

subplot(2, 1, 2);
plot(filtered_complex_signal);
title('Filtered Complex Signal');

% Point 6 (Frequency Response)
w_rect_wave = linspace(0, 2*pi, length(rect_wave));
w_complex_signal = linspace(0, 2*pi, length(complex_signal));

H_rect_wave = fft(rect_wave);
H_complex_signal = fft(real(complex_signal));

figure;
subplot(2, 2, 1);
plot(w_rect_wave, abs(H_rect_wave));
title('Magnitude Response - Rectangular Wave');
xlabel('Frequency (rad/sample)');
ylabel('Magnitude');

subplot(2, 2, 2);
plot(w_complex_signal, abs(H_complex_signal));
title('Magnitude Response - Complex Signal');
xlabel('Frequency (rad/sample)');
ylabel('Magnitude');

subplot(2, 2, 3);
plot(w_rect_wave, angle(H_rect_wave));
title('Phase Response - Rectangular Wave');
xlabel('Frequency (rad/sample)');
ylabel('Phase (rad)');

subplot(2, 2, 4);
plot(w_complex_signal, angle(H_complex_signal));
title('Phase Response - Complex Signal');
xlabel('Frequency (rad/sample)');
ylabel('Phase (rad)');

% 7 (Shifted DTFT)
% Define the rectangular wave
t_shifted = linspace(-pi, pi, 1000);
rect_wave = rectpuls(t_shifted, 1);

% Define a complex signal (for example, a sinusoidal signal)
complex_signal_shifted = exp(1j * 2 * pi * 3 * t_shifted); % Adjust the parameters as needed

% Compute DTFT of Rectangular Wave
samples_len = length(t_shifted);
omega = linspace(-pi, pi, samples_len);
DTFT_rect_wave = fftshift(fft(rect_wave));

% Plot the DTFT of Rectangular Wave
figure;
subplot(3, 2, 1);
plot(omega, abs(DTFT_rect_wave));
title('DTFT_Rectangular Wave')
xlabel('Omega');
ylabel('X*(e^{j\omega})');

% Compute DT FT of Complex Signal
DTFT_complex_signal_shifted = fftshift(fft(complex_signal_shifted));

% Plot the DTFT of Complex Signal
subplot(3, 2, 2);
plot(omega, abs(DTFT_complex_signal_shifted));
title('DTFT_Complex Signal')
xlabel('Omega');
ylabel('X*(e^{j\omega})');

% Shift the Time Vector and Compute Shifted DTFT for Rectangular Wave
omega_shift_rect = 60;
t_shifted_rect = linspace(-pi, pi, samples_len) - omega_shift_rect;
DTFT_rect_wave_shifted = fftshift(fft(rect_wave .* exp(-1j * omega_shift_rect * t_shifted_rect)));

% Plot the Shifted DTFT for Rectangular Wave
subplot(3, 2, 3);
plot(omega, abs(DTFT_rect_wave_shifted));
title('DTFT_Rectangular Wave_Shifted');
xlabel('Omega');
ylabel('X*(e^{j\omega})');

% Shift the Time Vector and Compute Shifted DTFT for Complex Signal
omega_shift_complex = 60; % Adjust the shift for the complex signal
t_shifted_complex = linspace(-pi, pi, samples_len) - omega_shift_complex;
DTFT_complex_signal_shifted = fftshift(fft(complex_signal_shifted .* exp(-1j * omega_shift_complex * t_shifted_complex)));

% Plot the Shifted DTFT for Complex Signal
subplot(3, 2, 4);
plot(omega, abs(DTFT_complex_signal_shifted));
title('DTFT_Complex Signal_Shifted');
xlabel('Omega');
ylabel('X*(e^{j\omega})');

% Convolution of DTFTs
convolution_result_DTFT = ifft(ifftshift(DTFT_rect_wave_shifted .* DTFT_complex_signal_shifted));

% Plot the Convolution of Shifted DTFTs
subplot(3, 2, 5);
plot(convolution_result_DTFT);
title('Convolution of Shifted DTFTs');
xlabel('Sample');
ylabel('Amplitude');

% Convolution in Time Domain for Comparison
convolution_result_time_domain = conv(rect_wave, complex_signal_shifted);

% Plot the Convolution in Time Domain for Comparison
subplot(3, 2, 6);
plot(convolution_result_time_domain(1:length(omega)));
title('Convolution in Time Domain (for Comparison)');
xlabel('Sample');
ylabel('Amplitude');

%---------------------------------------
% Point 8 (Discrete Fourier Transform (DFT) and Its Properties)
x_rect_wave = ones(size(t_shifted));
x_complex_signal = exp(1j * 2 * pi * freq * t_shifted);

figure;

% Compute and Plot DFT for Original Length
dft_result_rect_wave = fft(x_rect_wave);
dft_result_complex_signal = fft(x_complex_signal);

subplot(4, 4, 1);
plot(abs(dft_result_rect_wave));
title('DFT Result - Rectangular Wave (Original Length)');

subplot(4, 4, 2);
plot(abs(dft_result_complex_signal));
title('DFT Result - Complex Signal (Original Length)');

% Compute and Plot 4-point DFT
dft_4_point_rect_wave = fft(x_rect_wave, 4);
dft_4_point_complex_signal = fft(x_complex_signal, 4);

subplot(4, 4, 3);
plot(abs(dft_4_point_rect_wave));
title('4-point DFT - Rectangular Wave');

subplot(4, 4, 4);
plot(abs(dft_4_point_complex_signal));
title('4-point DFT - Complex Signal');

% Compute and Plot 8-point DFT
dft_8_point_rect_wave = fft(x_rect_wave, 8);
dft_8_point_complex_signal = fft(x_complex_signal, 8);

subplot(4, 4, 5);
plot(abs(dft_8_point_rect_wave));
title('8-point DFT - Rectangular Wave');

subplot(4, 4, 6);
plot(abs(dft_8_point_complex_signal));
title('8-point DFT - Complex Signal');

% Compute and Plot 16-point DFT
dft_16_point_rect_wave = fft(x_rect_wave, 16);
dft_16_point_complex_signal = fft(x_complex_signal, 16);

subplot(4, 4, 7);
plot(abs(dft_16_point_rect_wave));
title('16-point DFT - Rectangular Wave');

subplot(4, 4, 8);
plot(abs(dft_16_point_complex_signal));
title('16-point DFT - Complex Signal');

% Compute and Plot DFT with Zero Padding (e.g., zero-padding by 8 zeros)
zero_padded_length = length(t_shifted) + 8;
dft_zero_padded_rect_wave = fft(x_rect_wave, zero_padded_length);
dft_zero_padded_complex_signal = fft(x_complex_signal, zero_padded_length);

subplot(4, 4, 9);
plot(abs(dft_zero_padded_rect_wave));
title('DFT with Zero Padding - Rectangular Wave');

subplot(4, 4, 10);
plot(abs(dft_zero_padded_complex_signal));
title('DFT with Zero Padding - Complex Signal');

% Point 9 (Z-transform)
syms z n;
rect_wave_z_transform = ztrans(rect_wave, n, z); 
complex_signal_z_transform = ztrans(real(complex_signal), n, z);   

disp(rect_wave_z_transform);
disp(complex_signal_z_transform); 
disp("The code execution is complete.");
