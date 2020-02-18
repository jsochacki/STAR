function [GAIN] = set_Hammerstein_PA_output_power(input, output_power_dBm, pa_coefficients, GINIT, linear_distortion_coefficients)

Q = length(linear_distortion_coefficients);
half_length_Q = (Q + 1) / 2;
N = length(input);
Z0 = 50; %not necessary since calculations are relative but in case you want to know actual power
GAIN = GINIT;

P_TARGET = Z0 * 0.001 * power(10, output_power_dBm / 10);
%FIND PTARGET
lin_input_wld = cconv(input*GAIN, linear_distortion_coefficients);
lin_input_wld = lin_input_wld(1+(half_length_Q):end-(half_length_Q));
temp_NLIN = Memoryless_Polynomial_Amplifier(lin_input_wld, pa_coefficients);
TX_POWER_OUT_dBm_NLIN = 10*log10((temp_NLIN*temp_NLIN')/(N*Z0*0.001));
while P_TARGET > TX_POWER_OUT_dBm_NLIN
   GAIN = GAIN + 0.01;
   lin_input_wld = cconv(input*GAIN, linear_distortion_coefficients);
   lin_input_wld = lin_input_wld(1+(half_length_Q):end-(half_length_Q));
   temp_NLIN = Memoryless_Polynomial_Amplifier(lin_input_wld, pa_coefficients);
   TX_POWER_OUT_dBm_NLIN = 10*log10((temp_NLIN*temp_NLIN')/(N*Z0*0.001));
end
%FOUND PTARGET

end