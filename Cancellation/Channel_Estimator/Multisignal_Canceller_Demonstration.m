clear all
clc

%Simulation Settings
MODCOD = 1;
SYMBOLS_PER_SLOT = 200000;
equalizer_length = 257;
half_equalizer_length = (equalizer_length - 1) / 2;

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

%Add AWGN
tx_waveform_normalized = tx_waveform ./ sqrt(oversampling_rate * ((tx_waveform * tx_waveform') / length(tx_waveform)));
[n0, sigma] = generate_awgn_from_EsNo(tx_waveform_normalized, 0,  1000, oversampling_rate);
rx_signal = tx_waveform_normalized + n0;

channel_response = [zeros(1,100) 1 zeros(1, 11) 1 zeros(1,100)];
half_channel_response_length = (length(channel_response) - 1) / 2;
rx_signal_post_channel = cconv(rx_signal, channel_response);
rx_signal_post_channel = rx_signal_post_channel(1+(half_channel_response_length):end-(half_channel_response_length));

channel_response_estimate_coefficients = Least_Squares_Linear_Solution(rx_signal_post_channel, rx_signal, equalizer_length);
channel_response_estimate_coefficients = fliplr(channel_response_estimate_coefficients);

equalizer_coefficients = Least_Squares_Linear_Solution(rx_signal, rx_signal_post_channel, equalizer_length);
equalizer_coefficients = fliplr(equalizer_coefficients);

rx_signal_post_equalization = cconv(rx_signal_post_channel, equalizer_coefficients);
rx_signal_post_equalization = rx_signal_post_equalization(1+(half_equalizer_length):end-(half_equalizer_length));

%Receive Filtering
baseband_waveform = cconv(rx_signal_post_channel, fliplr(filter_h));
baseband_symbols = downsample(baseband_waveform(1+((2*ringing_length)):end-((2*ringing_length))), oversampling_rate);

SNR_dB_without_equalization = Measure_SNR(baseband_symbols, symbol_stream);
EVM_percent_without_equalization = 100*sqrt(1/power(10,SNR_dB_without_equalization/10));

%Receive Filtering
baseband_waveform_pd = cconv(rx_signal_post_equalization, fliplr(filter_h));
baseband_symbols_pd = downsample(baseband_waveform_pd(1+(2*ringing_length):end-(2*ringing_length)), oversampling_rate);

SNR_dB_with_equalization = Measure_SNR(baseband_symbols_pd, symbol_stream);
EVM_percent_with_equalization = 100*sqrt(1/power(10,SNR_dB_with_equalization/10));

%plot
[gain_figure gain_axis] = create_eq_gain_plot([], [], ...
                 10*log10(((abs(tx_waveform)).^2)/(length(tx_waveform)*50*0.001)), ...
                 10*log10(((abs(rx_signal_post_channel)).^2)/(length(rx_signal_post_channel)*50*0.001)) - 10*log10(((abs(tx_waveform)).^2)/(length(tx_waveform)*50*0.001)), ...
                 10*log10(((abs(tx_waveform)).^2)/(length(tx_waveform)*50*0.001)), ...
                 10*log10(((abs(rx_signal_post_equalization)).^2)/(length(rx_signal_post_equalization)*50*0.001)) - 10*log10(((abs(tx_waveform)).^2)/(length(tx_waveform)*50*0.001)), ...
                 [], ...
                 [], ...
                 -90, -30, -10, 10);

[psd_figure psd_axis] = create_eq_psd_plot([], [], ...
                linspace(-1/2, 1/2, length(tx_waveform)), ...
                20*log10(filtfilt(1./(500*ones(1,500)),1,abs(fftshift(fft(tx_waveform, length(tx_waveform)))))), ...
                linspace(-1/2, 1/2, length(rx_signal_post_channel)), ...
                20*log10(filtfilt(1./(500*ones(1,500)),1,abs(fftshift(fft(rx_signal_post_channel, length(rx_signal_post_channel)))))), ...
                linspace(-1/2, 1/2, length(rx_signal_post_equalization)), ...
                20*log10(filtfilt(1./(500*ones(1,500)),1,abs(fftshift(fft(rx_signal_post_equalization, length(rx_signal_post_equalization)))))), ...
                oversampling_rate, ...
                -30, ...
                70);

[constellation_figure constellation_axis] = create_eq_consteallation_plot([], [], ...
                           real(baseband_symbols), imag(baseband_symbols), ...
                           real(baseband_symbols_pd), imag(baseband_symbols_pd), ...
                           [], [], ...
                           SNR_dB_without_equalization, EVM_percent_without_equalization, ...
                           SNR_dB_with_equalization, EVM_percent_with_equalization, ...
                           -1.5, 1.5, -1.5, 1.5);
