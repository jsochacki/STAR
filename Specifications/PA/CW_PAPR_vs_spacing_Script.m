%%This is just to look at PAPR vs tone spacing
clear all
clc

%Simulation Settings
f1 = 1;
tone_spacing_percents = (0:0.1:100)/100;
oversampling_rate = 32;
SAMPLES = 31250 * oversampling_rate;
TONES = 3;

t = 0:1:(SAMPLES - 1);
t = t ./ oversampling_rate;

PAPRS = zeros(length(tone_spacing_percents), 1);
for n = 1:1:length(tone_spacing_percents)
   signal = sin(2*pi*f1*t);
   for nn = 2:1:TONES
      frequency = f1 * (1 + ((nn - 1) * tone_spacing_percents(n)));
      signal = signal + sin(2*pi*frequency*t);
   end
   PAPRS(n, 1) = PAPR_dB(signal, []);
end

plot(tone_spacing_percents, PAPRS)