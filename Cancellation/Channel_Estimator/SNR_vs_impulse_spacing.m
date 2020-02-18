clear all
clc

%Simulation Settings
MODCOD = 1;
SYMBOLS_PER_SLOT = 2000;
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

N_impulses = 2;
equalizer_lengths = [2.^(9:1:11)] + 1; 
for nn = 1:1:length(equalizer_lengths)
   equalizer_length = equalizer_lengths(nn);
   half_equalizer_length = (equalizer_length - 1) / 2;
   minimum_fb_padding = 2*(2*ringing_length);
   impulse_spacings{nn} = 1:2:(equalizer_length - (minimum_fb_padding + N_impulses));
   for n = 1:1:length(impulse_spacings{nn})
      pad = zeros(1,(equalizer_length - (impulse_spacings{nn}(n) + N_impulses)) / 2);
      channel_response = [pad 1 zeros(1, impulse_spacings{nn}(n)) power(10, -S_I_dB / 20) pad];
      length_channel_response = length(channel_response);
      half_channel_response_length = (length_channel_response - 1) / 2;
      rx_signal_post_channel = cconv(rx_signal, channel_response);
      rx_signal_post_channel = rx_signal_post_channel(1+(half_channel_response_length):end-(half_channel_response_length));

      %Receive Filtering
      baseband_waveform = cconv(rx_signal_post_channel, fliplr(filter_h));
      baseband_symbols = downsample(baseband_waveform(1+((2*ringing_length)):end-((2*ringing_length))), oversampling_rate);

      SNR_dB_without_equalization = Measure_SNR(baseband_symbols, symbol_stream);
      EVM_percent_without_equalization = 100*sqrt(1/power(10,SNR_dB_without_equalization/10));

      equalizer_coefficients = Least_Squares_Linear_Solution(rx_signal, rx_signal_post_channel, equalizer_length);
      equalizer_coefficients = fliplr(equalizer_coefficients);

      rx_signal_post_equalization = cconv(rx_signal_post_channel, equalizer_coefficients);
      rx_signal_post_equalization = rx_signal_post_equalization(1+(half_equalizer_length):end-(half_equalizer_length));

      %Receive Filtering
      baseband_waveform_pd = cconv(rx_signal_post_equalization, fliplr(filter_h));
      baseband_symbols_pd = downsample(baseband_waveform_pd(1+(2*ringing_length):end-(2*ringing_length)), oversampling_rate);

      length_ratio(nn) = length_channel_response/equalizer_length;
      SNR_dB_with_equalization(nn, n) = Measure_SNR(baseband_symbols_pd, symbol_stream);
      EVM_percent_with_equalization(nn, n) = 100*sqrt(1/power(10,SNR_dB_with_equalization(nn, n)/10));
   end
   length_ratio(nn) = length_channel_response/equalizer_length;
end

figure1 = figure;
axes1 = axes('Parent',figure1);
hold(axes1,'on');
ylabel('SNR (dB)');
xlabel('Echo Spacing in taps');
title(['Post Echo Cancellation SNR vs equalizer length', newline, 'with equalizer length to channel response ratio = 1']);
xlim(axes1,[min(impulse_spacings{end}) max(impulse_spacings{end})]);
ylim(axes1,[min(SNR_dB_with_equalization(end,:)) max(SNR_dB_with_equalization(end,:))]);
box(axes1,'on');
set(axes1,'XGrid','on','YGrid','on', 'XScale','log');
legend1 = legend(axes1,'show');
set(legend1,...
   'Position',[0.411317569996558 0.457936510898231 0.206081076223101 0.109730845899562]);
for nn = 1:1:length(equalizer_lengths)
   plot_vec = SNR_dB_with_equalization(nn,1:1:length(impulse_spacings{nn}));
   plot(axes1, impulse_spacings{nn}, plot_vec, ...
      'DisplayName',num2str([equalizer_lengths(nn)], 'equalizer length = %d'), 'Marker', '.');
      %'DisplayName',num2str([length_ratio(nn)], 'Equalizer length / channel response length ratio = %d'), 'Marker','o');
end
['Shannon Power Efficiency Limit',newline,'(Minimum Eb/No for Maximum Spectral Efficiency)']);