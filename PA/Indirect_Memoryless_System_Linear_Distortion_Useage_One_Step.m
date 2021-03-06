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
PAPR_Reduction = 1;

%Generate temporary PA coefficients
pa_coefficients = [14.9740 + j*0.519 0 -23.0954 + j*4.9680 0 21.3936 + j*0.4305 0];
hup = [0 0.1 0 0 0 1 0 0 0 0.1 0];

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

half_length_hup = (length(hup) - 1) / 2;
Q = filter_length_in_symbols * oversampling_rate + 1;
half_length_Q = (Q - 1) / 2;
linear_distortion_h = cconv(fir1((filter_length_in_symbols * oversampling_rate), 1 / oversampling_rate), hup);
linear_distortion_h = linear_distortion_h(1+(half_length_hup):1:(end-(half_length_hup)));
half_linear_distortion_length = (length(linear_distortion_h) - 1) / 2;
linear_distortion_length_in_symbols = (length(linear_distortion_h) - 1) / oversampling_rate;

shift = (find(max(filter_h)==filter_h)) - (find(max(linear_distortion_h)==linear_distortion_h));

%Generate Constellation and tx waveform
[Complex_Alphabet Binary_Alphabet Decimal_Alphabet BITS_PER_WORD] = dvbs2_Constellations(MODCOD);
symbol_stream = randsrc(1, SYMBOLS_PER_SLOT, Complex_Alphabet);
oversampled_symbol_stream = upsample(symbol_stream, oversampling_rate);
tx_waveform = cconv(oversampled_symbol_stream, filter_h);

%Set output back off and generate waveform at the output of the pa
tx_waveform_wld = cconv(tx_waveform, linear_distortion_h);
tx_waveform_wld = tx_waveform_wld((1+half_linear_distortion_length):1:(end-half_linear_distortion_length));
[REQUIRED_SIGNAL_GAIN SYSTEM_POWER_GAIN_dB] = set_OBO(tx_waveform, OBO_FROM_P1DB, pa_coefficients, 0.01);
%non observable but calculate here for information sake
[REQUIRED_SIGNAL_GAIN_WLD SYSTEM_POWER_GAIN_dB_WLD] = set_OBO(tx_waveform_wld, OBO_FROM_P1DB, pa_coefficients, 0.01);
tx_signal = tx_waveform*REQUIRED_SIGNAL_GAIN;
tx_signal_wld = tx_waveform_wld*REQUIRED_SIGNAL_GAIN;
tx_waveform_at_pa_output = Memoryless_Polynomial_Amplifier(tx_signal_wld, pa_coefficients);
%tx_waveform_at_pa_output = Memoryless_Polynomial_Amplifier(tx_signal*power(10, -PREDISTORTER_BACKOFF/10), pa_coefficients);
tx_power_at_pa_output_nopd = 10*log10((tx_waveform_at_pa_output*tx_waveform_at_pa_output')/(length(tx_waveform_at_pa_output)*50*0.001));

%Add AWGN
tx_waveform_at_pa_output_normalized = tx_waveform_at_pa_output ./ sqrt(oversampling_rate * ((tx_waveform_at_pa_output * tx_waveform_at_pa_output') / length(tx_waveform_at_pa_output)));
[n0, sigma] = generate_awgn_from_EsNo(tx_waveform_at_pa_output_normalized, 0,  1000, oversampling_rate);
rx_signal = tx_waveform_at_pa_output_normalized + n0;
%Receive Filtering
baseband_waveform = cconv(rx_signal, fliplr(filter_h));
baseband_symbols = downsample(baseband_waveform(1+((2*ringing_length)-shift):end-((2*ringing_length)+shift)), oversampling_rate);


SNR_dB_without_predistortion = Measure_SNR(baseband_symbols, symbol_stream);
EVM_percent_without_predistortion = 100*sqrt(1/power(10,SNR_dB_without_predistortion/10));

%Now Apply Ideal Predistortion
%linear estimated based on unobservable signals
%Take initial guess at b coefficients using subset of data so inversion is possible
[linear_distortion_coefficients X] = Least_Squares_Linear_Solution(tx_signal_wld, tx_signal, Q);
%Since this passes back the solution to the equation Y = X*b which is 
%convolution with fliplr(b) (i.e. correlation with b'.' since correlation
%should also have the complex conjugate taken)
%what b to be in the form of a filter for convolution not correlation
linear_distortion_coefficients = fliplr(linear_distortion_coefficients);
tx_waveform_wld_estimate = cconv(tx_waveform, linear_distortion_coefficients);
figure(5)
plot(abs(tx_waveform_wld((1 + half_length_Q):1:(end - (half_length_Q)))))
hold on
plot(abs(X.'*fliplr(linear_distortion_coefficients).'))
plot(abs(tx_waveform_wld_estimate((1 + (2*half_length_Q)):1:(end - (2*half_length_Q)))))

[equalizer_coefficients X] = Least_Squares_Linear_Solution(tx_waveform, tx_waveform_wld, Q);
tx_waveform_wold = cconv(tx_waveform_wld, equalizer_coefficients);
tx_waveform_wold = tx_waveform_wold((1+half_length_Q):1:(end-half_length_Q));

%linear estimated based on observable signals
[linear_distortion_coefficients X] = Least_Squares_Linear_Solution(tx_waveform_at_pa_output / power(10, SYSTEM_POWER_GAIN_dB/20), tx_signal, Q);
linear_distortion_coefficients = fliplr(linear_distortion_coefficients);
tx_waveform_wld_estimate = cconv(tx_signal, linear_distortion_coefficients);
figure(6)
plot(abs(tx_waveform_wld((1 + half_length_Q):1:(end - (half_length_Q)))))
hold on
plot(abs(X.'*fliplr(linear_distortion_coefficients).'))
plot(abs(tx_waveform_wld_estimate((1 + (2*half_length_Q)):1:(end - (2*half_length_Q)))))

%Now apply
[equalizer_coefficients X] = Least_Squares_Linear_Solution(tx_signal, tx_waveform_at_pa_output / power(10, SYSTEM_POWER_GAIN_dB/20), Q);

if CFR
   pre_CFR_PAPR = PAPR_dB(tx_signal, []);
   [tx_signal_pre_pd post_CFR_PAPR] = serial_peak_cancellation(tx_signal, filter_h, pre_CFR_PAPR - PAPR_Reduction, CFR_Iterations);
else
   tx_signal_pre_pd = tx_signal;
end

tx_waveform_equalized = cconv(tx_signal_pre_pd, equalizer_coefficients);
tx_waveform_equalized = tx_waveform_equalized((1+(half_length_Q-shift)):1:(end-(half_length_Q+shift)));
tx_waveform_eq_wld = cconv(tx_waveform_equalized, linear_distortion_h);
tx_waveform_eq_wld = tx_waveform_eq_wld((1+(half_linear_distortion_length-shift)):1:(end-(half_linear_distortion_length+shift)));

good_power = (1/length(tx_signal)) * sum(tx_signal .* conj(tx_signal));
[tx_waveform_eq_wld AGC_GAIN] = AGC(tx_waveform_eq_wld, good_power);

tx_waveform_at_pa_output_pre_pd = Memoryless_Polynomial_Amplifier(tx_waveform_eq_wld, pa_coefficients);

pd_coefficients = Least_Squares_Memoryless_Odd_Polynomial_Solution(tx_waveform_eq_wld, tx_waveform_at_pa_output_pre_pd / power(10, SYSTEM_POWER_GAIN_dB/20), POLYNOMIAL_ORDER);
pd_tx_waveform = Memoryless_Polynomial_Amplifier(tx_waveform_eq_wld*power(10, -PREDISTORTER_BACKOFF/20), pd_coefficients);

tx_waveform_at_pa_output_pd = Memoryless_Polynomial_Amplifier(pd_tx_waveform, pa_coefficients);

tx_power_at_pa_output_pd = 10*log10((tx_waveform_at_pa_output_pd*tx_waveform_at_pa_output_pd')/(length(tx_waveform_at_pa_output_pd)*50*0.001));

%Add AWGN
%Have to normailze to symbol power here before adding noise so you don't
%ruin the mssp off the rx signal in the SNR measurement
tx_waveform_at_pa_output_pd_normalized = tx_waveform_at_pa_output_pd ./ sqrt(oversampling_rate * ((tx_waveform_at_pa_output_pd * tx_waveform_at_pa_output_pd') / length(tx_waveform_at_pa_output_pd)));
[n0, sigma] = generate_awgn_from_EsNo(tx_waveform_at_pa_output_pd_normalized, 0,  1000, oversampling_rate);
rx_signal_pd = tx_waveform_at_pa_output_pd_normalized + n0;

%Receive Filtering
baseband_waveform_pd = cconv(rx_signal_pd, fliplr(filter_h));
baseband_symbols_pd = downsample(baseband_waveform_pd(1+((2*ringing_length)-shift):end-((2*ringing_length)+shift)), oversampling_rate);


SNR_dB_with_predistortion = Measure_SNR(baseband_symbols_pd, symbol_stream);
EVM_percent_with_predistortion = 100*sqrt(1/power(10,SNR_dB_with_predistortion/10));

%plot
[gain_figure gain_axis] = create_gain_plot([], [], ...
                 10*log10(((abs(tx_signal)).^2)/(length(tx_signal)*50*0.001)), ...
                 10*log10(((abs(tx_waveform_at_pa_output)).^2)/(length(tx_waveform_at_pa_output)*50*0.001)) - 10*log10(((abs(tx_signal)).^2)/(length(tx_signal)*50*0.001)), ...
                 10*log10(((abs(tx_signal)).^2)/(length(tx_signal)*50*0.001)), ...
                 10*log10(((abs(tx_waveform_at_pa_output_pd)).^2)/(length(tx_waveform_at_pa_output_pd)*50*0.001)) - 10*log10(((abs(tx_signal)).^2)/(length(tx_signal)*50*0.001)), ...
                 [10*log10(((min(abs(tx_signal)).^2))/(length(tx_signal)*50*0.001)) 10*log10(((max(abs(tx_signal)).^2))/(length(tx_signal)*50*0.001))], ...
                 SYSTEM_POWER_GAIN_dB*[1 1], ...
                 -75, -50, 23, 31);

[psd_figure psd_axis] = create_psd_plot([], [], ...
                linspace(-1/2, 1/2, length(tx_signal)), ...
                20*log10(filtfilt(1./(100*ones(1,100)),1,abs(fftshift(fft(tx_signal, length(tx_signal)))))), ...
                linspace(-1/2, 1/2, length(tx_waveform_at_pa_output)), ...
                20*log10(filtfilt(1./(100*ones(1,100)),1,abs(fftshift(fft(tx_waveform_at_pa_output, length(tx_waveform_at_pa_output)))))), ...
                linspace(-1/2, 1/2, length(tx_waveform_at_pa_output_pd)), ...
                20*log10(filtfilt(1./(100*ones(1,100)),1,abs(fftshift(fft(tx_waveform_at_pa_output_pd, length(tx_waveform_at_pa_output_pd)))))), ...
                oversampling_rate, ...
                10, ...
                80);

[constellation_figure constellation_axis] = create_consteallation_plot([], [], ...
                           real(baseband_symbols), imag(baseband_symbols), ...
                           real(baseband_symbols_pd), imag(baseband_symbols_pd), ...
                           real(symbol_stream), imag(symbol_stream), ...
                           SNR_dB_without_predistortion, EVM_percent_without_predistortion, ...
                           SNR_dB_with_predistortion, EVM_percent_with_predistortion, ...
                           -1.5, 1.5, -1.5, 1.5);