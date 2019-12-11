clear all
clc

%Simulation Settings
MODCOD = 1; %QPSK
% MODCOD = 13; %32 APSK
% MODCOD = 19; %DVBT2 16 QAM
SYMBOLS_PER_SLOT = 200000;

%generate tx and rx filters
oversampling_rate = 8;
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

tx_waveform_power = 10*log10((tx_waveform*tx_waveform')/(length(tx_waveform)*50*0.001));

%Add AWGN
tx_waveform_normalized = tx_waveform ./ sqrt(oversampling_rate * ((tx_waveform * tx_waveform') / length(tx_waveform)));
[n0, sigma] = generate_awgn_from_EsNo(tx_waveform_normalized, 0,  1000, oversampling_rate);
tx_waveform_normalized_with_AWGN = tx_waveform_normalized + n0;

tx_signal = tx_waveform_normalized_with_AWGN;
pre_CFR_PAPR = PAPR_dB(tx_signal, []);
SNR_dB_without_predistortion = Measure_SNR(symbol_stream, symbol_stream);
EVM_percent_without_predistortion = 100*sqrt(1/power(10,SNR_dB_without_predistortion/10));
SNR_dB_with_predistortion = -100;
EVM_percent_with_predistortion = 0;

nWrites = write_aeroflex_file( tx_signal, './QPSK_4x_oversampling_6p25_alpha_200ksymbols_NO_AWGN_7dB_PAPR.aiq', false );

%plot
[gain_figure gain_axis] = create_gain_plot([], [], ...
                 10*log10(((abs(tx_signal)).^2)/(length(tx_signal)*50*0.001)), ...
                 10*log10(((abs(tx_signal)).^2)/(length(tx_signal)*50*0.001)) - 10*log10(((abs(tx_signal)).^2)/(length(tx_signal)*50*0.001)), ...
                 [], ...
                 [], ...
                 [], ...
                 [], ...
                 -75, -50, 23, 31);

[psd_figure psd_axis] = create_psd_plot([], [], ...
                linspace(-1/2, 1/2, length(tx_signal)), ...
                20*log10(filtfilt(1./(100*ones(1,100)),1,abs(fftshift(fft(tx_signal, length(tx_signal)))))), ...
                [], ...
                [], ...
                [], ...
                [], ...
                oversampling_rate, ...
                10, ...
                80);

[constellation_figure constellation_axis] = create_consteallation_plot([], [], ...
                           real(symbol_stream), imag(symbol_stream), ...
                           [], [], ...
                           [], [], ...
                           SNR_dB_without_predistortion, EVM_percent_without_predistortion, ...
                           SNR_dB_with_predistortion, EVM_percent_with_predistortion, ...
                           -1.5*max(abs(symbol_stream)), 1.5*max(abs(symbol_stream)), -1.5*max(abs(symbol_stream)), 1.5*max(abs(symbol_stream)));
