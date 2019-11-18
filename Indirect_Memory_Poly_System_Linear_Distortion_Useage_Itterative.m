clear all
clc

%Simulation Settings
MODCOD = 1;
SYMBOLS_PER_SLOT = 200000;
OBO_FROM_P1DB = 0;
PREDISTORTER_BACKOFF = 1;
CFR = 0;
CFR_Iterations = 100;
QPA = 3; %must be odd
PA_POLYNOMIAL_ORDER = 5;
QPD = 5; %must be odd
DPD_POLYNOMIAL_ORDER = 7;

%Generate temporary PA coefficients
pa_coefficients = [ 1.4991 + j*0.0173  4.9913 + j*0.1730  0.0991 + j*0.0173 ...
                   -0.0698 - j*0.1656 -7.6984 - j*1.6566 -0.7698 - j*0.1656 ...
                    0.7131 - j*0.0143  7.1312 - j*0.1435  0.0131 - j*0.0143] * 4;

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

half_length_QPA = (QPA - 1) / 2;

%Generate Constellation and tx waveform
[Complex_Alphabet Binary_Alphabet Decimal_Alphabet BITS_PER_WORD] = dvbs2_Constellations(MODCOD);
symbol_stream = randsrc(1, SYMBOLS_PER_SLOT, Complex_Alphabet);
oversampled_symbol_stream = upsample(symbol_stream, oversampling_rate);
tx_waveform = cconv(oversampled_symbol_stream, filter_h);
base_signal_PAPR_dB = PAPR_dB(tx_waveform, []);

%Set output back off and generate waveform at the output of the pa
[REQUIRED_SIGNAL_GAIN SYSTEM_POWER_GAIN_dB] = set_memory_PA_OBO(tx_waveform, OBO_FROM_P1DB, pa_coefficients, 0.01, QPA);
tx_signal = tx_waveform*REQUIRED_SIGNAL_GAIN;
tx_waveform_at_pa_output = Memory_Polynomial_Amplifier(tx_signal, pa_coefficients, PA_POLYNOMIAL_ORDER, QPA);
%tx_waveform_at_pa_output = Memoryless_Polynomial_Amplifier(tx_signal*power(10, -PREDISTORTER_BACKOFF/10), pa_coefficients);
tx_power_at_pa_output_nopd = 10*log10((tx_waveform_at_pa_output*tx_waveform_at_pa_output')/(length(tx_waveform_at_pa_output)*50*0.001));

%get pa coefficients with memory
pa_coefficients_with_memory = Memory_Polynomial_Solution(tx_waveform_at_pa_output, tx_signal, PA_POLYNOMIAL_ORDER, QPA);
tx_waveform_at_pa_output_modeled = Memory_Polynomial_Amplifier(tx_signal, pa_coefficients_with_memory, PA_POLYNOMIAL_ORDER, QPA);

figure(1)
plot([10*log10(((min(abs(tx_signal)).^2))/(length(tx_signal)*50*0.001)) 10*log10(((max(abs(tx_signal)).^2))/(length(tx_signal)*50*0.001))], ...
     SYSTEM_POWER_GAIN_dB+[10*log10(((min(abs(tx_signal)).^2))/(length(tx_signal)*50*0.001)) 10*log10(((max(abs(tx_signal)).^2))/(length(tx_signal)*50*0.001))], 'k')
hold on, grid on
plot(10*log10(((abs(tx_signal)).^2)/(length(tx_signal)*50*0.001)), ...
     10*log10(((abs(tx_waveform_at_pa_output)).^2)/(length(tx_signal)*50*0.001)), 'b.')
plot(10*log10(((abs(tx_signal)).^2)/(length(tx_signal)*50*0.001)), ...
     10*log10(((abs(tx_waveform_at_pa_output_modeled)).^2)/(length(tx_signal)*50*0.001)), 'k.')

figure(2)
plot([10*log10(((min(abs(tx_signal)).^2))/(length(tx_signal)*50*0.001)) 10*log10(((max(abs(tx_signal)).^2))/(length(tx_signal)*50*0.001))], ...
     SYSTEM_POWER_GAIN_dB*[1 1], 'k')
hold on, grid on
plot(10*log10(((abs(tx_signal)).^2)/(length(tx_signal)*50*0.001)), ...
     10*log10(((abs(tx_waveform_at_pa_output)).^2)/(length(tx_signal)*50*0.001)) - 10*log10(((abs(tx_signal)).^2)/(length(tx_signal)*50*0.001)), 'b.')
plot(10*log10(((abs(tx_signal)).^2)/(length(tx_signal)*50*0.001)), ...
     10*log10(((abs(tx_waveform_at_pa_output_modeled)).^2)/(length(tx_signal)*50*0.001)) - 10*log10(((abs(tx_signal)).^2)/(length(tx_signal)*50*0.001)), 'k.')
  
figure(3)
plot(20*log10(filtfilt(1./(100*ones(1,100)),1,abs(fftshift(fft(tx_signal, length(tx_signal)))))))
hold on
plot(20*log10(filtfilt(1./(100*ones(1,100)),1,abs(fftshift(fft(tx_waveform_at_pa_output, length(tx_signal)))))))
plot(20*log10(filtfilt(1./(100*ones(1,100)),1,abs(fftshift(fft(tx_waveform_at_pa_output_modeled, length(tx_waveform_at_pa_output_modeled)))))))

%Add AWGN
tx_waveform_at_pa_output_normalized = tx_waveform_at_pa_output ./ sqrt(oversampling_rate * ((tx_waveform_at_pa_output * tx_waveform_at_pa_output') / length(tx_waveform_at_pa_output)));
[n0, sigma] = generate_awgn_from_EsNo(tx_waveform_at_pa_output_normalized, 0,  1000, oversampling_rate);
rx_signal = tx_waveform_at_pa_output_normalized + n0;

%Receive Filtering
baseband_waveform = cconv(rx_signal, fliplr(filter_h));
baseband_symbols = downsample(baseband_waveform(1+((2*ringing_length)):end-((2*ringing_length))), oversampling_rate);
figure(4)
plot(baseband_symbols, 'bo')
hold on, grid on

SNR_dB_without_predistortion = Measure_SNR(baseband_symbols, symbol_stream);
EVM_percent_without_predistortion = 100*sqrt(1/power(10,SNR_dB_without_predistortion/10));

%Now solve for a single itteration of the memory prolynomial predistorter
training_signal = tx_signal;
tx_waveform_at_pa_output_pd = tx_waveform_at_pa_output;
for n = 1:1:10
   %Get initial memoryless pd coefficients
   pd_coefficients = Memory_Polynomial_Solution(training_signal, tx_waveform_at_pa_output_pd / power(10, SYSTEM_POWER_GAIN_dB/20), DPD_POLYNOMIAL_ORDER, QPD);
   pd_tx_waveform = Memory_Polynomial_Amplifier(tx_signal*power(10, -PREDISTORTER_BACKOFF/20), pd_coefficients, DPD_POLYNOMIAL_ORDER, QPD);
   
   if CFR
      pre_CFR_PAPR = PAPR_dB(pd_tx_waveform, []);
      [pd_tx_waveform_post_cfr post_CFR_PAPR] = serial_peak_cancellation(pd_tx_waveform, filter_h, base_signal_PAPR_dB, CFR_Iterations);
      tx_waveform_at_pa_output_pd = Memory_Polynomial_Amplifier(pd_tx_waveform_post_cfr, pa_coefficients, PA_POLYNOMIAL_ORDER, QPA);
      training_signal = pd_tx_waveform_post_cfr;
   else
      tx_waveform_at_pa_output_pd = Memory_Polynomial_Amplifier(pd_tx_waveform, pa_coefficients, PA_POLYNOMIAL_ORDER, QPA);
      training_signal = pd_tx_waveform;
   end
   
   if n == 9 PREDISTORTER_BACKOFF = 1.8; end
   tx_power_at_pa_output_pd = 10*log10((tx_waveform_at_pa_output_pd*tx_waveform_at_pa_output_pd')/(length(tx_waveform_at_pa_output_pd)*50*0.001));
end

figure(1)
plot(10*log10(((abs(tx_signal).^2)/(length(tx_signal)*50*0.001))), ...
     10*log10(((abs(tx_waveform_at_pa_output_pd)).^2)/(length(tx_signal)*50*0.001)), 'r.')

figure(2)
plot(10*log10(((abs(tx_signal)).^2)/(length(tx_signal)*50*0.001)), ...
     10*log10(((abs(tx_waveform_at_pa_output_pd)).^2)/(length(tx_signal)*50*0.001)) - 10*log10(((abs(tx_signal)).^2)/(length(tx_signal)*50*0.001)), 'r.')
  
figure(3)
plot(20*log10(filtfilt(1./(100*ones(1,100)),1,abs(fftshift(fft(tx_waveform_at_pa_output_pd, length(tx_signal)))))))

%Add AWGN
%Have to normailze to symbol power here before adding noise so you don't
%ruin the mssp off the rx signal in the SNR measurement
tx_waveform_at_pa_output_pd_normalized = tx_waveform_at_pa_output_pd ./ sqrt(oversampling_rate * ((tx_waveform_at_pa_output_pd * tx_waveform_at_pa_output_pd') / length(tx_waveform_at_pa_output_pd)));
[n0, sigma] = generate_awgn_from_EsNo(tx_waveform_at_pa_output_pd_normalized, 0,  1000, oversampling_rate);
rx_signal_pd = tx_waveform_at_pa_output_pd_normalized + n0;

%Receive Filtering
baseband_waveform_pd = cconv(rx_signal_pd, fliplr(filter_h));
baseband_symbols_pd = downsample(baseband_waveform_pd(1+(2*ringing_length):end-(2*ringing_length)), oversampling_rate);
figure(4)
plot(baseband_symbols_pd, 'ro')

SNR_dB_with_predistortion = Measure_SNR(baseband_symbols_pd, symbol_stream);
EVM_percent_with_predistortion = 100*sqrt(1/power(10,SNR_dB_with_predistortion/10));