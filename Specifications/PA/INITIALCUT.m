clear all
clc

%Simulation Settings
MODCOD = 1;
SYMBOLS_PER_SLOT = 200000;
OBO_FROM_P1DB = 0;
PREDISTORTER_BACKOFF = 1.8;
CFR = 0;
CFR_Iterations = 100;
Waveform_PAPR = 1;
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

if CFR
   pre_CFR_PAPR = PAPR_dB(tx_waveform, []);
   [tx_waveform post_CFR_PAPR] = serial_peak_cancellation(tx_waveform, filter_h, Waveform_PAPR, CFR_Iterations);
   tx_waveform = tx_waveform;
else
   tx_waveform = tx_waveform;
end

%Set output back off and generate waveform at the output of the pa
[REQUIRED_SIGNAL_GAIN SYSTEM_POWER_GAIN_dB] = set_memory_PA_OBO(tx_waveform, OBO_FROM_P1DB, pa_coefficients, 0.01, QPA);
tx_signal = tx_waveform*REQUIRED_SIGNAL_GAIN;
tx_waveform_at_pa_output = Memory_Polynomial_Amplifier(tx_signal, pa_coefficients, PA_POLYNOMIAL_ORDER, QPA);
tx_power_at_pa_output = 10*log10((tx_waveform_at_pa_output*tx_waveform_at_pa_output')/(length(tx_waveform_at_pa_output)*50*0.001));


%Add AWGN
tx_waveform_at_pa_output_normalized = tx_waveform_at_pa_output ./ sqrt(oversampling_rate * ((tx_waveform_at_pa_output * tx_waveform_at_pa_output') / length(tx_waveform_at_pa_output)));
[n0, sigma] = generate_awgn_from_EsNo(tx_waveform_at_pa_output_normalized, 0,  1000, oversampling_rate);
rx_signal = tx_waveform_at_pa_output_normalized + n0;

%Receive Filtering
baseband_waveform = cconv(rx_signal, fliplr(filter_h));
baseband_symbols = downsample(baseband_waveform(1+((2*ringing_length)):end-((2*ringing_length))), oversampling_rate);


SNR_dB = Measure_SNR(baseband_symbols, symbol_stream);
EVM_percent = 100*sqrt(1/power(10,SNR_dB/10));

%plot
[gain_figure gain_axis] = create_gain_plot([], [], ...
                 10*log10(((abs(tx_signal)).^2)/(length(tx_signal)*50*0.001)), ...
                 10*log10(((abs(tx_signal)).^2)/(length(tx_signal)*50*0.001)) - 10*log10(((abs(tx_signal)).^2)/(length(tx_signal)*50*0.001)), ...
                 10*log10(((abs(tx_signal)).^2)/(length(tx_signal)*50*0.001)), ...
                 10*log10(((abs(tx_waveform_at_pa_output)).^2)/(length(tx_waveform_at_pa_output)*50*0.001)) - 10*log10(((abs(tx_signal)).^2)/(length(tx_signal)*50*0.001)), ...
                 [], ...
                 [], ...
                 -75, -50, 23, 31);

[psd_figure psd_axis] = create_psd_plot([], [], ...
                linspace(-1/2, 1/2, length(tx_signal)), ...
                20*log10(filtfilt(1./(100*ones(1,100)),1,abs(fftshift(fft(tx_signal, length(tx_signal)))))), ...
                linspace(-1/2, 1/2, length(tx_waveform_at_pa_output)), ...
                20*log10(filtfilt(1./(100*ones(1,100)),1,abs(fftshift(fft(tx_waveform_at_pa_output, length(tx_waveform_at_pa_output)))))), ...
                [], ...
                [], ...
                oversampling_rate, ...
                10, ...
                80);

[constellation_figure constellation_axis] = create_consteallation_plot([], [], ...
                           real(baseband_symbols), imag(baseband_symbols), ...
                           [], [], ...
                           real(Complex_Alphabet), imag(Complex_Alphabet), ...
                           SNR_dB, EVM_percent, ...
                           [], [], ...
                           -1.5, 1.5, -1.5, 1.5);
