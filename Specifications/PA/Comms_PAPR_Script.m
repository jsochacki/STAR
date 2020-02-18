clear all
clc

%Simulation Settings
MODCODS = [1:1:14 19];
alphas = [0.3 0.25 0.2 0.15 0.01 0.05];
SYMBOLS_PER_SLOT = 1000000;

oversampling_rate = 4;
filter_length_in_symbols = 48;
filter_implementation_type = 'firrcoswu';  
filter_half_filter_length_at_design_rate = (filter_length_in_symbols .* oversampling_rate) / 2;
ringing_length = filter_half_filter_length_at_design_rate;

PAPRS = zeros(length(MODCODS), length(alphas));
for n = 1:1:length(alphas)
   for nn = 1:1:length(MODCODS)
      %generate tx and rx filters
      
      [filter_h, result] = generate_srrc_filter(filter_implementation_type, ...
                                                filter_length_in_symbols, ...
                                                alphas(n), ...
                                                oversampling_rate);

      %Makes unity gain filter
      filter_h = filter_h ./ sqrt(sum(power(filter_h, 2)));

      %Generate Constellation and tx waveform
      [Complex_Alphabet] = dvbs2_Constellations(MODCODS(nn));
      symbol_stream = randsrc(1, SYMBOLS_PER_SLOT, Complex_Alphabet);
      oversampled_symbol_stream = upsample(symbol_stream, oversampling_rate);
      tx_waveform = cconv(oversampled_symbol_stream, filter_h);
      PAPRS(nn, n) = PAPR_dB(tx_waveform, []);
   end
end