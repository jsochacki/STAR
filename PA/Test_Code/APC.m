function [output_after_apc POWER_GAIN_dB] = APC(signal, desired_power_dBm)
%msp = @(x) (1/length(x)) * sum(x .* x'.');
%A = msp(signal);

P_dBm = (1/(50*0.001)) * (1/length(signal)) * sum(signal .* conj(signal));
output_after_apc = signal * sqrt(desired_power_dBm / P_dBm);
POWER_GAIN_dB = 10*log10(desired_power_dBm / P_dBm);

end