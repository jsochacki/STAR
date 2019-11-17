clear all
clc

%Simulation Settings
MODCOD = 1;
SYMBOLS_PER_SLOT = 200000;
OBO_FROM_P1DB = 0;
POLYNOMIAL_ORDER = 5;
PREDISTORTER_BACKOFF = 1;
CFR = 0;
CFR_Iterations = 100;

%Generate temporary PA coefficients
x = [1 2 3 4 5 6 7 8 9 10];
y = [2 4 6 7.5 9 10 11 11.3 11.7 12] + j*[0.001 0.01 0.03 0.06 0.1 0.2 0.3 0.4 0.4 0.3];
pa_coefficients = Least_Squares_Memoryless_Odd_Polynomial_Solution(y, x, POLYNOMIAL_ORDER);
pa_coefficients = [14.9740 + j*0.519 0 -23.0954 + j*4.9680 0 21.3936 + j*0.4305 0];

%generate tx and rx filters
oversampling_rate = 4;
filter_alpha = 0.0625;
filter_length_in_symbols = 48;
filter_implementation_type = 'firrcoswu';

filter_half_filter_length_at_design_rate = (filter_length_in_symbols .* oversampling_rate) / 2;
ringing_length = filter_half_filter_length_at_design_rate;

[filter_h, result] = generate_srrc_filter(filter_implementation_type, ...
                                          filter_length_in_symbols, ...
                                          filter_alpha, ...
                                          oversampling_rate);

%Makes unity gain filter
filter_h = filter_h ./ sqrt(sum(power(filter_h, 2)));

%Generate Constellation and tx waveform
[Complex_Alphabet Binary_Alphabet Decimal_Alphabet BITS_PER_WORD] = dvbs2_Constellations(MODCOD);
symbol_stream = randsrc(1, SYMBOLS_PER_SLOT, Complex_Alphabet);
oversampled_symbol_stream = upsample(symbol_stream, oversampling_rate);
tx_waveform = cconv(oversampled_symbol_stream, filter_h);
base_signal_PAPR_dB = PAPR_dB(tx_waveform, []);

%Set output back off and generate waveform at the output of the pa
[REQUIRED_SIGNAL_GAIN SYSTEM_POWER_GAIN_dB] = set_OBO(tx_waveform, OBO_FROM_P1DB, pa_coefficients, 0.01);
tx_signal = tx_waveform*REQUIRED_SIGNAL_GAIN;
tx_waveform_at_pa_output = Memoryless_Polynomial_Amplifier(tx_signal, pa_coefficients);
%tx_waveform_at_pa_output = Memoryless_Polynomial_Amplifier(tx_signal*power(10, -PREDISTORTER_BACKOFF/10), pa_coefficients);
tx_power_at_pa_output_nopd = 10*log10((tx_waveform_at_pa_output*tx_waveform_at_pa_output')/(length(tx_waveform_at_pa_output)*50*0.001));

figure(1)
plot([10*log10(((min(abs(tx_signal)).^2))/(length(tx_signal)*50*0.001)) 10*log10(((max(abs(tx_signal)).^2))/(length(tx_signal)*50*0.001))], ...
     SYSTEM_POWER_GAIN_dB+[10*log10(((min(abs(tx_signal)).^2))/(length(tx_signal)*50*0.001)) 10*log10(((max(abs(tx_signal)).^2))/(length(tx_signal)*50*0.001))], 'k')
hold on, grid on
plot(10*log10(((abs(tx_signal)).^2)/(length(tx_signal)*50*0.001)), ...
     10*log10(((abs(tx_waveform_at_pa_output)).^2)/(length(tx_waveform_at_pa_output)*50*0.001)), 'b.')

figure(2)
plot([10*log10(((min(abs(tx_signal)).^2))/(length(tx_signal)*50*0.001)) 10*log10(((max(abs(tx_signal)).^2))/(length(tx_signal)*50*0.001))], ...
     SYSTEM_POWER_GAIN_dB*[1 1], 'k')
hold on, grid on
plot(10*log10(((abs(tx_signal)).^2)/(length(tx_signal)*50*0.001)), ...
     10*log10(((abs(tx_waveform_at_pa_output)).^2)/(length(tx_waveform_at_pa_output)*50*0.001)) - 10*log10(((abs(tx_signal)).^2)/(length(tx_signal)*50*0.001)), 'b.')
  
figure(3)
plot(20*log10(filtfilt(1./(100*ones(1,100)),1,abs(fftshift(fft(tx_signal, length(tx_signal)))))))
hold on
plot(20*log10(filtfilt(1./(100*ones(1,100)),1,abs(fftshift(fft(tx_waveform_at_pa_output, length(tx_waveform_at_pa_output)))))))

%Add AWGN
[n0, sigma] = generate_awgn_from_EsNo(tx_waveform_at_pa_output, 0,  1000, oversampling_rate);
rx_signal = tx_waveform_at_pa_output + n0;

%Receive Filtering
baseband_waveform = cconv(rx_signal, fliplr(filter_h));
baseband_symbols = downsample(baseband_waveform, oversampling_rate);
baseband_symbols = baseband_symbols(1+(filter_length_in_symbols):end-(filter_length_in_symbols));
figure(4)
plot(baseband_symbols, 'bo')
hold on, grid on


%Now Apply Ideal Predistortion
pd_coefficients = Least_Squares_Memoryless_Odd_Polynomial_Solution(tx_signal, tx_waveform_at_pa_output / power(10, SYSTEM_POWER_GAIN_dB/20), POLYNOMIAL_ORDER);
%pd_coefficients = Least_Squares_Memoryless_Odd_Polynomial_Solution(tx_signal, tx_waveform_at_pa_output, POLYNOMIAL_ORDER);
pd_tx_waveform = Memoryless_Polynomial_Amplifier(tx_signal*power(10, -PREDISTORTER_BACKOFF/20), pd_coefficients);
%pd_tx_waveform = Memoryless_Polynomial_Amplifier(tx_signal, pd_coefficients);
if CFR
   pre_CFR_PAPR = PAPR_dB(pd_tx_waveform, []);
   [pd_tx_waveform_post_cfr post_CFR_PAPR] = serial_peak_cancellation(pd_tx_waveform, filter_h, base_signal_PAPR_dB, CFR_Iterations);
   tx_waveform_at_pa_output_pd = Memoryless_Polynomial_Amplifier(pd_tx_waveform_post_cfr, pa_coefficients);
else
   tx_waveform_at_pa_output_pd = Memoryless_Polynomial_Amplifier(pd_tx_waveform, pa_coefficients);
end
tx_power_at_pa_output_pd = 10*log10((tx_waveform_at_pa_output_pd*tx_waveform_at_pa_output_pd')/(length(tx_waveform_at_pa_output_pd)*50*0.001));
figure(1)
plot(10*log10(((abs(tx_signal)).^2)/(length(tx_signal)*50*0.001)), ...
     10*log10(((abs(tx_waveform_at_pa_output_pd)).^2)/(length(tx_waveform_at_pa_output_pd)*50*0.001)), 'r.')

figure(2)
plot(10*log10(((abs(tx_signal)).^2)/(length(tx_signal)*50*0.001)), ...
     10*log10(((abs(tx_waveform_at_pa_output_pd)).^2)/(length(tx_waveform_at_pa_output_pd)*50*0.001)) - 10*log10(((abs(tx_signal)).^2)/(length(tx_signal)*50*0.001)), 'r.')
  
figure(3)
plot(20*log10(filtfilt(1./(100*ones(1,100)),1,abs(fftshift(fft(tx_waveform_at_pa_output_pd, length(tx_waveform_at_pa_output_pd)))))))

%Add AWGN
[n0, sigma] = generate_awgn_from_EsNo(tx_waveform_at_pa_output_pd, 0,  1000, oversampling_rate);
rx_signal_pd = tx_waveform_at_pa_output_pd + n0;

%Receive Filtering
baseband_waveform_pd = cconv(rx_signal_pd, fliplr(filter_h));
baseband_symbols_pd = downsample(baseband_waveform_pd, oversampling_rate);
baseband_symbols_pd = baseband_symbols_pd(1+(filter_length_in_symbols):end-(filter_length_in_symbols));
figure(4)
plot(baseband_symbols_pd, 'ro')

%Now Apply Ideal Predistortion
pd_coefficients2 = Least_Squares_Memoryless_Odd_Polynomial_Solution(tx_signal*power(10, -PREDISTORTER_BACKOFF/20), tx_waveform_at_pa_output_pd / power(10, SYSTEM_POWER_GAIN_dB/20), POLYNOMIAL_ORDER);
%pd_coefficients = Least_Squares_Memoryless_Odd_Polynomial_Solution(tx_signal, tx_waveform_at_pa_output, POLYNOMIAL_ORDER);
pd_tx_waveform = Memoryless_Polynomial_Amplifier(tx_signal*power(10, -PREDISTORTER_BACKOFF/20), pd_coefficients2);
pd_tx_waveform2 = Memoryless_Polynomial_Amplifier(pd_tx_waveform, pd_coefficients);
pd_coefficients_cumulative = Least_Squares_Memoryless_Odd_Polynomial_Solution(tx_signal*power(10, -PREDISTORTER_BACKOFF/20), pd_tx_waveform2 / power(10, 0/20), POLYNOMIAL_ORDER);
pd_tx_waveform2f = Memoryless_Polynomial_Amplifier(tx_signal*power(10, -PREDISTORTER_BACKOFF/20), pd_coefficients_cumulative);
%pd_tx_waveform = Memoryless_Polynomial_Amplifier(tx_signal, pd_coefficients);
if CFR
   pre_CFR_PAPR = PAPR_dB(pd_tx_waveform2f, []);
   [pd_tx_waveform_post_cfr post_CFR_PAPR] = serial_peak_cancellation(pd_tx_waveform2f, filter_h, base_signal_PAPR_dB, CFR_Iterations);
   tx_waveform_at_pa_output_pd2 = Memoryless_Polynomial_Amplifier(pd_tx_waveform_post_cfr, pa_coefficients);
else
   tx_waveform_at_pa_output_pd2 = Memoryless_Polynomial_Amplifier(pd_tx_waveform2f, pa_coefficients);
end
tx_power_at_pa_output_pd = 10*log10((tx_waveform_at_pa_output_pd2*tx_waveform_at_pa_output_pd2')/(length(tx_waveform_at_pa_output_pd2)*50*0.001));
figure(1)
plot(10*log10(((abs(tx_signal)).^2)/(length(tx_signal)*50*0.001)), ...
     10*log10(((abs(tx_waveform_at_pa_output_pd2)).^2)/(length(tx_waveform_at_pa_output_pd2)*50*0.001)), 'c.')

figure(2)
plot(10*log10(((abs(tx_signal)).^2)/(length(tx_signal)*50*0.001)), ...
     10*log10(((abs(tx_waveform_at_pa_output_pd2)).^2)/(length(tx_waveform_at_pa_output_pd2)*50*0.001)) - 10*log10(((abs(tx_signal)).^2)/(length(tx_signal)*50*0.001)), 'c.')
  
figure(3)
plot(20*log10(filtfilt(1./(100*ones(1,100)),1,abs(fftshift(fft(tx_waveform_at_pa_output_pd2, length(tx_waveform_at_pa_output_pd2)))))))

%Add AWGN
[n0, sigma] = generate_awgn_from_EsNo(tx_waveform_at_pa_output_pd2, 0,  1000, oversampling_rate);
rx_signal_pd = tx_waveform_at_pa_output_pd2 + n0;

%Receive Filtering
baseband_waveform_pd = cconv(rx_signal_pd, fliplr(filter_h));
baseband_symbols_pd = downsample(baseband_waveform_pd, oversampling_rate);
baseband_symbols_pd = baseband_symbols_pd(1+(filter_length_in_symbols):end-(filter_length_in_symbols));
figure(4)
plot(baseband_symbols_pd, 'co')
