clear all
clc

%Simulation Settings
MODCOD = 1;
SYMBOLS_PER_SLOT = 2000;
equalizer_lengths = [2.^(6:1:11)] + 1; 
S_I_dB = 10;

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

channel_response = [zeros(1,100) 1 zeros(1, 11) power(10, -S_I_dB / 20) zeros(1,100)];
length_channel_response = length(channel_response);
half_channel_response_length = (length_channel_response - 1) / 2;
rx_signal_post_channel = cconv(rx_signal, channel_response);
rx_signal_post_channel = rx_signal_post_channel(1+(half_channel_response_length):end-(half_channel_response_length));

%Receive Filtering
baseband_waveform = cconv(rx_signal_post_channel, fliplr(filter_h));
baseband_symbols = downsample(baseband_waveform(1+((2*ringing_length)):end-((2*ringing_length))), oversampling_rate);

SNR_dB_without_equalization = Measure_SNR(baseband_symbols, symbol_stream);
EVM_percent_without_equalization = 100*sqrt(1/power(10,SNR_dB_without_equalization/10));

for n = 1:1:length(equalizer_lengths)
   half_equalizer_length = (equalizer_lengths(n) - 1) / 2;

   equalizer_coefficients = Least_Squares_Linear_Solution(rx_signal, rx_signal_post_channel, equalizer_lengths(n));
   equalizer_coefficients = fliplr(equalizer_coefficients);

   rx_signal_post_equalization = cconv(rx_signal_post_channel, equalizer_coefficients);
   rx_signal_post_equalization = rx_signal_post_equalization(1+(half_equalizer_length):end-(half_equalizer_length));

   %Receive Filtering
   baseband_waveform_pd = cconv(rx_signal_post_equalization, fliplr(filter_h));
   baseband_symbols_pd = downsample(baseband_waveform_pd(1+(2*ringing_length):end-(2*ringing_length)), oversampling_rate);

   SNR_dB_with_equalization(n) = Measure_SNR(baseband_symbols_pd, symbol_stream);
   EVM_percent_with_equalization(n) = 100*sqrt(1/power(10,SNR_dB_with_equalization(n)/10));
end

figure1 = figure;
axes1 = axes('Parent',figure1);
hold(axes1,'on');
ylabel('SNR (dB)');
xlabel('Equalizer Taps');
title('Post Echo Cancellation SNR vs equalizer length');
xlim(axes1,[min(equalizer_lengths) max(equalizer_lengths)]);
ylim(axes1,[min(SNR_dB_with_equalization) max(SNR_dB_with_equalization)]);
box(axes1,'on');
set(axes1,'XGrid','on','YGrid','on');
plot(axes1, equalizer_lengths, SNR_dB_with_equalization,'Marker','o');

figure2 = figure;
axes1 = axes('Parent',figure2);
hold(axes1,'on');
ylabel('SNR (dB)');
xlabel('Equalizer Taps / channel response length');
title('Post Echo Cancellation SNR vs equalizer length to channel response ratio');
xlim(axes1,[min(equalizer_lengths/length_channel_response) max(equalizer_lengths/length_channel_response)]);
ylim(axes1,[min(SNR_dB_with_equalization) max(SNR_dB_with_equalization)]);
box(axes1,'on');
set(axes1,'XGrid','on','YGrid','on');
plot(axes1, equalizer_lengths/length_channel_response, SNR_dB_with_equalization,'Marker','o');
