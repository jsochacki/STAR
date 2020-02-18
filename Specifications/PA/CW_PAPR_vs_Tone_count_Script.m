clear all
clc

%This is to look at PAPR vs number of tones for a single spacing (it will
%be the same for ALL spacings)

%Simulation Settings
f1 = 1;
tone_spacing_percent = 0.32;
oversampling_rate = 32;
SAMPLES = 31250 * oversampling_rate;
TONES = 40;

t = 0:1:(SAMPLES - 1);
t = t ./ oversampling_rate;

PAPRS = zeros(TONES, 1);
for n = 1:1:TONES
   signal = sin(2*pi*f1*t);
   for nn = 2:1:n
      frequency = f1 * (1 + ((nn - 1) * tone_spacing_percent));
      signal = signal + sin(2*pi*frequency*t);
   end
   PAPRS(n, 1) = PAPR_dB(signal, []);
end

plot(1:1:TONES, PAPRS, 'b-')
hold on
plot(1:1:TONES, PAPRS, 'r.')
